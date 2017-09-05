# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:03:58 2016

@author: Tobias Jachowski
"""
import persistent
from collections import deque


class Node(persistent.Persistent):
    """
    State information via references:

    Node implements a directed acyclic graph with depth and breadth first
    traversal. Each node has an attribute `cargo`, which can hold any python
    object. Node provides methods to get and set nodes (and cargo) based on
    different conditions, like level of depth or type of instance.
    """

    def __init__(self, cargo=None, **kwargs):
        """
        Attributes
        ----------
        cargo : object
            The cargo of the Node.
        """
        # initialize variables
        self._relatives = {'children': [],
                           'parents': []}
        self._max_relatives = {'children': -1,
                               'parents': -1}
        self.cargo = cargo
        # self._processed = None
        super().__init__(**kwargs)

    def add_child(self, child, index=None, after=None, before=None):
        """
        Add a child. See `add_relative()` for details.
        """
        return self.add_relative(child, child=True, index=index, after=after,
                                 before=before)

    def add_parent(self, parent, index=None, after=None, before=None):
        """
        Add a parent. See `add_relative()` for details.
        """
        return self.add_relative(parent, child=False, index=index, after=after,
                                 before=before)

    def add_relative(self, relative, child=True, index=None, after=None,
                     before=None, bidirectional=True):
        """
        Add a `relative` as a child or a parent.

        The index of the position of inserting the relative is determined by
        the following parameters in descending superseding order: `index`,
        `after`, and `before`.

        Parameters
        ----------
        relative : Node or cargo
            Node to be added as a relative, or cargo with attribute `_node`
            referencing a Node.
        child : bool, optional
            Add a child or a parent.
        index : int, optional
        after : Node or cargo, optional
        before : Node or cargo, optional
        bidirectional : bool, optional
            Add alsoe `self` to the `relative`, after `relative` has been added
            to `self`, i.e. not only reference the `relative` from `self`, but
            also reference `self` from `relative`.

        Returns
        -------
        bool
            True if `relative` could be added or `relative` is already present.
            Otherwise return False.
        """
        # Get nodes of relative, after, and before
        relative = _get_node(relative)
        if after is not None:
            after = _get_node(after)
        if before is not None:
            before = _get_node(before)

        # Try to add relative
        self_add = self._add_relative(relative, child=child, index=index,
                                      after=after, before=before)

        # Try to add self to relative
        relation_add = True
        if bidirectional:
            relation_add = relative._add_relative(self, child=not child)

        # check for circular reference
        circular_reference = relative.circular_reference(descendants=child)

        # Relative could not be added, or self could not be added to relative,
        # or circular_reference detected
        if not self_add or not relation_add or circular_reference:
            # Either relative couldn't be added or
            # circular reference occured -> remove relative
            self._remove_relative(relative, child=child)
            relative._remove_relative(self, child=not child)
            return False

        # relative was added or already present
        return True

    def _add_relative(self, relative, child=True, index=None, after=None,
                      before=None):
        if child:
            relatives = self._children
            max_relatives = self.max_children
        else:
            relatives = self._parents
            max_relatives = self.max_parents

        if max_relatives != -1 and len(relatives) >= max_relatives:
            print("Can not add relative. No more relatives allowed!")
            return False

        if relative is not None:
            if relative in relatives:
                # parent already exists, avoid adding a second time and return
                # True
                print("Relative %s already referenced." % relative)
                return True

            # determine the index of self._parents/self._childs before which
            # the parent/child should be inserted
            if index is not None:
                # highest priority of index determination, take the given index
                pass
            elif after in relatives:
                # second highest priority (determine index by after)
                index = relatives.index(after) + 1
            elif before in relatives:
                # third highest priority (determine index by before)
                index = relatives.index(before)
            else:
                # default: append parent/child (add parent/child at the end)
                index = len(relatives)

            # Add parent/child at determined index
            relatives.insert(index, relative)
            # inform ZODB of change
            self._p_changed = True

            return True

    def remove_child(self, child):
        return self.remove_relative(child, child=True)

    def remove_parent(self, parent):
        return self.remove_relative(parent, child=False)

    def remove_relative(self, relative, child=True):
        """
        Remove child or parent from self.
        """
        relative = _get_node(relative)

        self_remove = self._remove_relative(relative, child)
        relate_remove = relative._remove_relative(self, child=not child)

        # No relative was removed
        if not (self_remove or relate_remove):
            return False

        # relative was removed
        return True

    def _remove_relative(self, relative, child=True):
        if child:
            relatives = self._children
        else:
            relatives = self._parents

        # remove child/parent from relatives
        if relative in relatives:
            relatives.remove(relative)
            # inform ZODB of change
            self._p_changed = True
            return True

        return False

    def circular_reference(self, descendants=True):
        """
        Check for circular references.
        """
        # determine direction of search for circular reference
        relatives = self.relatives(descendants=descendants)

        # Search for a circular reference.
        # First, assume there is no circular reference:
        # If this View has no modifications or relatives, there is no
        # circular reference and circular stays False
        # Check for circular reference via modifications
        # for modification in modifications:
        #     circular = circular or modification.circular_reference(down,
        #                                                            caller)
        # Check for circular reference via relatives
        for relative in relatives:
            if relative is self:
                print("Circular reference detected -> Illegal connection!")
                return True

        # stop recursion: return result of circular recursive search
        return False

    def relatives(self, descendants=True, includeself=False, dft=True,
                  level=-1, cargo=False):
        """
        Traverse relatives.

        Parameters
        ----------
        descendants : bool
            Yield descendants or ancestors
        includeself : bool
            Yield only relatives or also include self
        dft : bool
            Use depth first or breadth first traversal
        level : int
            Yield self (0), up to first generation (1), up to second generation
            (2), ...
        cargo : bool
            Yield node or its cargo

        Yields
        -------
        Node
            If `cargo` is False.
        object
            If `cargo` is True.

        Notes
        -----
        Level describes the generation up to which a relative should be
        yielded:
                         |
        level 0         Node
                         |  \
        level 1          |   Node
                         |    |
        level 2          |   Node
                         |  / |  \
        level 1/3       Node Node Node
                       / |  \
        level 2/4  Node Node Node
                              |
        """
        # see http://jeremykun.com/2013/01/22/depth-and-breadth-first-search/

        # Initialization ...
        toprocess = deque()

        if dft:
            # reversed iteration for depth first search
            iterswitch = reversed
            addtoprocess = toprocess.append
        else:
            iterswitch = iter
            addtoprocess = toprocess.appendleft

        # start with either self or the relatives
        if includeself:
            # add self at level 0
            toprocess.append((0, self))
        elif level == -1 or level >= 1:
            # relatives are either children or parents:
            if descendants:
                relatives = self._children
            else:
                relatives = self._parents

            for n in iterswitch(relatives):
                # add relatives at level 1
                addtoprocess((1, n))

        visited = set()
        # better execution time but not thread save!
        # search_time = datetime.now().timestamp()
        # process_id = uuid.uuid4()
        # search_time = process_id
        # better execution time but not thread save!

        # Traversal and processing ...
        while len(toprocess) > 0:
            current_level, node = toprocess.pop()

            # if node._processed != search_time:
            if node not in visited:
                visited.add(node)
                # node._processed = search_time

                if cargo and node.cargo is not None:
                    yield node.cargo
                if not cargo:
                    yield node

                if descendants:
                    relatives = node._children
                else:
                    relatives = node._parents

                next_level = current_level + 1
                if level == -1 or next_level <= level:

                    for n in iterswitch(relatives):
                        # if n._processed != search_time:
                        if n not in visited:
                            addtoprocess((next_level, n))

    @property
    def _children(self):
        return self._relatives['children']

    @property
    def _parents(self):
        return self._relatives['parents']

    @property
    def max_children(self):
        return self._max_relatives['children']

    @max_children.setter
    def max_children(self, max_children):
        self._max_relatives['children'] = max_children

    @property
    def max_parents(self):
        return self._max_relatives['parents']

    @max_parents.setter
    def max_parents(self, max_parents):
        self._max_relatives['parents'] = max_parents

    @property
    def parents(self):
        return self.relatives(descendants=False, dft=False, level=1)

    @property
    def children(self):
        return self.relatives(descendants=True, dft=False, level=1)

    @property
    def ancestors(self):
        return self.relatives(descendants=False, dft=False)

    @property
    def descendants(self):
        return self.relatives(descendants=True, dft=False)


class GraphMember(persistent.Persistent):
    """
    GraphMember has one instance of class Node, offering convenience functions
    to handle the Node simply by inheriting, and give the ancestor class a Node
    like behaviour.
    It implements methods to inform other members of the graph of a change.
    """
    def __init__(self, max_parents=-1, max_children=-1, updated=False,
                 name=None, group=None, **kwargs):
        self._node = Node(cargo=self)

        self.max_parents = max_parents
        self.max_children = max_children

        # Initialize the status of a GraphMember to be not updated, per default
        self.updated = updated

        self.name = name
        self.group = group

    def members(self, name=None, group=None, instance_class=None,
                descendants=True, includeself=True, dft=False, level=-1):
        for relative in self._node.relatives(descendants=descendants,
                                             includeself=includeself, dft=dft,
                                             level=level, cargo=True):
            if (name is None or relative.name == name) \
                    and (group is None or relative.group == group) \
                    and (instance_class is None
                         or isinstance(relative, instance_class)):
                yield relative

    def group_ends(self, group, root=None, descendants=True, includeself=True):
        for member in self.members(group=group, descendants=descendants,
                                   includeself=includeself, dft=True):
            if member.is_end(root=root):
                yield member

    def group_roots(self, group, descendants=True, includeself=True):
        return self.group_ends(group, root=True, descendants=descendants,
                               includeself=includeself)

    def group_leafs(self, group, descendants=True, includeself=True):
        return self.group_ends(group, root=False, descendants=descendants,
                               includeself=includeself)

    def group_end(self, group, root=True, descendants=True, includeself=True):
        ends = self.group_ends(group, root=root, descendants=descendants,
                               includeself=includeself)
        return next(ends, None)

    def group_root(self, group, descendants=True, includeself=True):
        return self.group_end(group, root=True, descendants=descendants,
                              includeself=includeself)

    def group_leaf(self, group, descendants=True, includeself=True):
        return self.group_end(group, root=False, descendants=descendants,
                              includeself=includeself)

    def is_end(self, root=None):
        """
        Parameters
        ----------
        root : None or bool, optional
            If root is None, check whether GraphMember is root or leaf of
            own group. If root is True, check whether Graphmember is root of
            group, if root is False check Graphmember beeing a leaf of group.
        """
        if root is None:
            return self.is_group_root or self.is_group_leaf
        elif root:
            return self.is_group_root
        else:
            return self.is_group_leaf

    @property
    def is_group_root(self):
        parent = self.parent
        return parent is None or self.group != parent.group

    @property
    def is_group_leaf(self):
        child = self.child
        return child is None or self.group != child.group

    def set_changed(self, descendants=True, includeself=True, level=-1,
                    **kwargs):
        """
        Inform all descendants/ancestors (Nodes that have this Node as a
        parent/child) about a change, so that the children inform their
        children, about the change.

        Has to be called upon any change of `self` and parents.

        It calls `member_changed(ancestor=descendants, **kwargs)` on `self`,
        if `includeself` is True, and on all descendants.

        Parameters
        ----------
        descendants : bool
            Inform descendants or ancestors about change
        level : int
            Up to which generation the change should be proclaimed
        includeself : bool
            Should self be informed, too, i.e. should member_changed() be
            called?
        **kwargs
            See member_changed() for parameters
        """
        # inform node about change
        # the node in turn will call member_changed() method of the cargos
        # requested (see set_changed() method of node)
        if includeself:
            self.member_changed(ancestor=descendants, calledfromself=True,
                                **kwargs)

        # Get either descendants or ancestors to be informed of the change
        members = self.members(descendants=descendants, includeself=False,
                               dft=False, level=level)

        for member in members:
            member.member_changed(ancestor=descendants, calledfromself=False,
                                  **kwargs)

    def member_changed(self, ancestor=True, calledfromself=False,
                       updated=False, **kwargs):
        # `self` triggered the change. Set updated status of `self` according
        # to the parameter `updated`.
        if calledfromself:
            self.updated = updated
        # An ancestor triggered the change and `self` is set to be outdated. A
        # change of descendants will be ignored.
        if not calledfromself and ancestor:
            self.updated = False

    def add_member(self, member, child=True, index=None, after=None,
                   before=None, set_changed=True, **kwargs):
        # Try to add the new member
        added = self._node.add_relative(member, child=child, index=index,
                                        after=after, before=before, **kwargs)
        # Inform the child or self and the children about change
        if added and set_changed:
            if child:
                # Inform the child about the addition of a new parent (i.e. the
                # addition of self)
                member.set_changed()
            else:
                # inform self and children about addition of a new parent
                self.set_changed()
        return added

    def add_parent(self, parent, index=None, after=None, before=None,
                   set_changed=True, **kwargs):
        # Try to add the new parent
        return self.add_member(parent, child=False, index=index, after=after,
                               before=before, set_changed=set_changed,
                               **kwargs)

    def add_child(self, child, index=None, after=None, before=None,
                  set_changed=True, **kwargs):
        # Try to add the new child
        return self.add_member(child, child=True, index=index, after=after,
                               before=before, set_changed=set_changed,
                               **kwargs)

    def remove_member(self, member, child=True, set_changed=True):
        # Try to remove the member
        removed = self._node.remove_relative(member, child=child)

        # Inform the child or self and the children about change
        if removed and set_changed:
            if child:
                # Inform the child about the loss of a parent (i.e. the loss
                # of self)
                member.set_changed()
            else:
                # Inform self and children about loss of a parent
                self.set_changed()

        return removed

    def remove_parent(self, parent, set_changed=True):
        """
        Remove parent from self.parents
        """
        return self.remove_member(parent, child=False, set_changed=set_changed)

    def remove_child(self, child, set_changed=True):
        """
        Remove child from self.children
        """
        return self.remove_member(child, child=True, set_changed=set_changed)

    def set_member(self, member, child=True, set_changed=True):
        """
        Replace all children or parents with `member`.
        """
        if child:
            members = self.children
        else:
            members = self.parents

        # Remove all members
        for old_member in members:
            self.remove_member(old_member, child=child, set_changed=set_changed)

        return self.add_member(member, child=child, set_changed=set_changed)

    def set_parent(self, parent, set_changed=True):
        """
        Replace all parents with `parent`.
        """
        return self.set_member(parent, child=False, set_changed=set_changed)

    def set_child(self, child, set_changed=True):
        """
        Replace all children with `child`.
        """
        return self.set_member(child, child=True, set_changed=set_changed)

    def first_ancestor_instance(self, instance_class, dft=True, level=-1):
        ancestors = self.members(instance_class=instance_class,
                                 descendants=False, includeself=False, dft=dft,
                                 level=level)
        return next(ancestors, None)

    def parent_instances(self, instance_class):
        members = self.members(instance_class=instance_class,
                               descendants=False, includeself=False, level=1)
        return members

    def child_instances(self, instance_class=None):
        members = self.members(instance_class=instance_class,
                               descendants=True, includeself=False, level=1)
        return members

    @property
    def parents(self):
        return self.members(descendants=False, includeself=False, level=1)

    @property
    def children(self):
        return self.members(descendants=True, includeself=False, level=1)

    @property
    def parent(self):
        """
        Return first parent or None.
        """
        return next(self.parents, None)

    @property
    def child(self):
        """
        Return first child or None.
        """
        return next(self.children, None)

    @property
    def max_parents(self):
        return self._node.max_parents

    @max_parents.setter
    def max_parents(self, max_parents):
        self._node.max_parents = max_parents

    @property
    def max_children(self):
        return self._node.max_children

    @max_children.setter
    def max_children(self, max_children):
        self._node.max_children = max_children

    def __str__(self):
        if self.name is None or self.group is None:
            return self.__repr__()
        else:
            return "".join(["".join(self.name).ljust(22),
                            "".join(self.group).ljust(18),
                            "".join(["<", self.__class__.__name__,
                                     ">"]).ljust(14)
                            ])


def _get_node(cargo):
    node = cargo
    if not isinstance(cargo, Node):
        try:
            node = cargo._node
        except:
            raise AttributeError("The cargo does not have an attribute `_node`"
                                 " of type Node, which is necessary to be used"
                                 " as a cargo from an instance of class Node.")
    return node
