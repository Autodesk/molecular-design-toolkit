import sys
from .exports import exports


@exports
def can_use_widgets():
    """ Expanded from from http://stackoverflow.com/a/34092072/1958900
    """
    if 'IPython' not in sys.modules:
        # IPython hasn't been imported, definitely not
        return False
    from IPython import get_ipython

    # check for `kernel` attribute on the IPython instance
    if getattr(get_ipython(), 'kernel', None) is None:
        return False

    try:
        import ipywidgets as ipy
        import traitlets
    except ImportError:
        return False

    if int(ipy.__version__.split('.')[0]) < 6:
        print 'WARNING: widgets require ipywidgets 6.0 or later'
        return False

    return True


LAYOUT_PROPS = set(("height width max_height max_width min_height min_width "
                    "visibility display overflow overflow_x overflow_y  border margin "
                    "padding top left bottom right order flex_flow align_items "
                    "flex align_self align_content justify_content").split())


@exports
def process_widget_kwargs(kwargs):
    from ipywidgets import Layout

    layout = kwargs.get('layout', None)

    for arg in kwargs.keys():
        if arg in LAYOUT_PROPS:
            if not layout:
                layout = Layout()
            setattr(layout, arg, kwargs.pop(arg))

    if layout is not None:
        kwargs['layout'] = layout

    return kwargs
