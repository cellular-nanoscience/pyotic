# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:06:44 2016

@author: Tobias Jachowski
"""
import collections
import cloudpickle
import inspect
import persistent
from ipywidgets import widgets
from IPython.display import display  # Used to display widgets in the notebook


class unboundfunction(persistent.Persistent):
    """
    Class to hold references to functions and still be able to pickle them.
    To reference the function you want to bind:
    function_reference = boundmethod(function)
    """
    def __init__(self, ft, **kwargs):
        self.ft = ft

    def __getstate__(self):
        return cloudpickle.dumps(self.ft)

    def __setstate__(self, cloudepickled_ft):
        self.ft = cloudpickle.loads(cloudepickled_ft)

    def __call__(self, *a, **kw):
        return self.ft(*a, **kw)


class boundmethod(persistent.Persistent):
    """
    Class to hold references to methods and still be able to pickle them.
    To reference the method you want to bind:
    method_reference = boundmethod(method)
    """
    def __init__(self, mt, **kwargs):
        self.mt = mt

    def __getstate__(self):
        return self.mt.__self__, self.mt.__func__.__name__

    def __setstate__(self, tuple_self_function_name):
        (s, fn) = tuple_self_function_name
        self.mt = getattr(s, fn)

    def __call__(self, *a, **kw):
        return self.mt(*a, **kw)


class Attributes(persistent.Persistent):
    """
    Store a dictionary in ZODB and use keys of the dict as attributes of an
    instance of Attributes.
    """
    def __init__(self):
        self._store = {}

    def values(self):
        return self._store.values()

    def __len__(self):
        return self._store.__len__()

    def __iter__(self):
        return self._store.__iter__()

    def __contains__(self, key):
        return key in self._store

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __setitem__(self, key, value):
        self.__setattr__(key, value)

    def __delitem__(self, key):
        del self._store[key]

    def __getattr__(self, name):
        """
        Allow keys of status to be used as attributes
        """
        if name in self._store:
            return self._store[name]
        else:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        A __setattr__() hook will also override the Persistent __setattr__()
        hook.  User code must treat it much like __getattribute__().  The
        user-defined code must call _p_setattr() first to handle special
        attributes; _p_setattr() takes the attribute name and value.  If it
        returns True, Persistent handled the attribute. If not, the user code
        can run. If the user code modifies the object’s state, it must assigned
        to _p_changed
        """
        if not self._p_setattr(name, value):
            if '_store' in self.__dict__:
                self._store[name] = value
                self._p_changed = True  # inform ZODB of change
            else:
                super().__setattr__(name, value)


class FakeWidget(object):
    """
    FakeWidget is used for environments, where the graphical backend does not
    support instantiation of ipywidgets.widgets (self._create_widget()).

    Parameters
    ----------
    value : int, float, bool or None
      The value of the FakeWidget. If the value is set to None, a FakeWidget
      whose value is set to None and whose value cannot be changed anymore is
      created, to represent a "FakeButton" widget.
    description : str
      This value should describe what the value is supposed to represent.
    """
    def __init__(self, value, description, **kwargs):
        self._value = value
        self.description = description

    def __repr__(self):
        value = self.value
        if value is None:
            value = 'Button'
        return ''.join((self.description, ': %s' % value))

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        # prevent changing the value from/to None
        if self._value is not None and value is not None:
            self._value = value

# make relation of value type and widget class
WIDGET_CLASS = {
    None: widgets.Button,
    bool: widgets.Checkbox,
    float: widgets.FloatText,
    int: widgets.FloatText,
    list: widgets.SelectMultiple,
    tuple: widgets.SelectMultiple
}


class InteractiveAttributes(persistent.Persistent):
    """
    Stores widgets and implements methods to add widgets, change values of
    widgets, access values of widgets, and link widgets to callback functions
    called upon value change.
    """

    def __init__(self, **kwargs):
        # Widgets which graphically represent values and thereby add
        # interactive functionality to change and inspect values
        self._widgets = {}
        # self._widgets = AttributeDict()

        # List of callback functions [function1, function2, ...], which are
        # called whenever a value of a widget is changed
        self._callback_functions = {}

        # Controls whether callback functions upon change of widget value are
        # called or not.
        # Set to False, when you don't want the callback functions called.
        self._widget_callbacks_active = True

    def display(self):
        """
        Display all widgets, usually in an ipython notebook.
        """
        for widget in self._widgets.values():
            display(widget)

    def add(self, key, value=None, description=None, callback_functions=None,
            options=None, **kwargs):
        """
        Adds a parameter, which can be represented as a widget.
        Key is the name of the widget.
        The widget is added to the modification and a callback function is
        registered which is called upon any change of the widget.

        callback_functions: array (default [])
            array containing functions (callback functions) to be called upon
            change of the parameter.
        """
        # create and register a widget
        widget = self._create_widget(key, value=value, description=description,
                                     options=options)
        self._widgets[key] = widget

        # Register callback functions to be called upon widget value change
        self._callback_functions[key] = callback_functions or []

        # replace callback_functions by picklable versions
        for i, function in enumerate(self._callback_functions[key]):
            self._callback_functions[key][i] = boundmethod(function)

        # inform ZODB of change
        self._p_changed = True

    def _create_widget(self, key, value=None, description=None, options=None):
        # create widget according to the type of value
        try:
            if value is None:
                value_type = None
            else:
                value_type = type(value)
            if isinstance(value, collections.Iterable) \
                and not isinstance(options, collections.Iterable):
                    options = value
            widget_class = WIDGET_CLASS[value_type]
            widget = widget_class(description=description, value=value,
                                  options=options)

            # register _widget_callback() to be called upon change of widget
            # value
            value_changed = self.__make_widget_callback(self._widget_callback,
                                                        key)
            if value is None:
                widget.on_click(value_changed)
            else:
                widget.observe(value_changed, names='value')
        except:
            # Create FakeWidget in case a graphical widget could not be created
            widget = FakeWidget(value, description)

        return widget

    def __make_widget_callback(self, function, key):
        # needed, see keyword "late binding" (here for key) or
        # http://stackoverflow.com/questions/3431676/creating-functions-in-a-loop
        # Called, whenever one of the values of the widgets is changed.
        def widget_value_changed(name, value=None):
            self._p_changed = True  # inform ZODB of change
            if self._widget_callbacks_active:
                function(key)
        return widget_value_changed

    def _widget_callback(self, key):
        """
        Set a param with key key in self._values to value value.
        call_triggers can be set to False to prevent subsequent calls of
        """
        self._call_cb_functions(key)

    def _call_cb_functions(self, key, **kwargs):
        # call callback_functions
        for function in self._callback_functions[key]:
            if inspect.getargspec(function.mt)[2] is None:
                function()
            else:
                function(**kwargs)

    def set_value(self, key, value, callback=True, **kwargs):
        """
        Sets the value of the widget with key 'key'. Callback functions
        (callback_functions) registered to the widget can be disabled by
        setting callback to False.
        """
        if key in self and value != self._widgets[key].value:
            # disable automatic widget callback
            self._widget_callbacks_active = False

            # set value of widget
            self._widgets[key].value = value
            self._p_changed = True  # inform ZODB of change
            # Is also called by widget_value_changed() in
            # __make_widget_callback() and could be neglected, because setting
            # the immutable value self._widget_callbacks_active already
            # triggers a self._p_changed = True

            # restore automatic widget callback
            self._widget_callbacks_active = True

            # call _callback_functions (what _widget_callback() would do) if
            # callback is wanted
            if callback:
                self._call_cb_functions(key, **kwargs)

    def __contains__(self, key):
        return key in self._widgets

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __setitem__(self, key, value):
        self.__setattr__(key, value)

    def __getattr__(self, name):
        """
        Allow attributes to be used as params selections for set/get_params
        """
        if name in self:
            return self._widgets[name].value
        else:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        A __setattr__() hook will also override the Persistent __setattr__()
        hook.  User code must treat it much like __getattribute__().  The
        user-defined code must call _p_setattr() first to handle special
        attributes; _p_setattr() takes the attribute name and value.  If it
        returns True, Persistent handled the attribute. If not, the user code
        can run. If the user code modifies the object’s state, it must assigned
        to _p_changed
        """
        if not self._p_setattr(name, value):
            if '_widgets' in self.__dict__ and name in self:
                self.set_value(name, value)
                # is done by set_value:
                # self._p_changed = True # inform ZODB of change
            else:
                super().__setattr__(name, value)

    def __getstate__(self):
        """Return state values to be pickled."""
        # To circumvent: TypeError: can't pickle instancemethod objects
        # __getstate__ should replace instances in self._widgets with
        # information necessary to recreate widgets:
        # widget_classe, description, and value
        state = self.__dict__.copy()
        state['_widgets'] = state['_widgets'].copy()
        for key, widget in state['_widgets'].items():
            value = None
            if hasattr(widget, 'value'):
                value = widget.value
            description = widget.description
            options = None
            if hasattr(widget, 'options'):
                options = widget.options
            # since 0.6.3 options are stored. For compatibility reasons and
            # future optional parameters, additionally store None as a 4th one
            state['_widgets'][key] = [value, description, options, None]
        return state

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        # __setstate__ should refill self._widgets with the information stored
        # in self._widgets. Compare with self.add() for necessary steps
        for key, widgetState in state['_widgets'].items():
            if len(widgetState) == 2:
                value, description = widgetState
                options = None
            elif len(widgetState) == 3:  # pre 0.5.0
                widget_class, description, value = widgetState
                if widget_class is widgets.Button:
                    value = None
                options = None
            else:  # len == 4, since 0.6.3
                value, description, options, _ = widgetState
            widget = self._create_widget(key, value=value, options=options,
                                         description=description)
            state['_widgets'][key] = widget
        self.__dict__.update(state)
