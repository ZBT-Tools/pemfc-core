from . import widget_set
from . import button
from . import frame


class WidgetFactory:

    def __init__(self):
        self._widget_set_factory = widget_set.WidgetSetFactory()
        self._button_factory = button.ButtonFactory()
        self._frame_factory = frame.FrameFactory()

    def create(self, frame, **kwargs):
        widget_type = kwargs.get('type', None)
        if widget_type is None:
            if 'sub_frame_dicts' in kwargs or 'widget_dicts' in kwargs:
                return self._frame_factory.create_frame(frame, **kwargs)
        else:
            if 'Set' in widget_type:
                return self._widget_set_factory.create(frame, **kwargs)
            elif 'Button' in widget_type and 'Set' not in widget_type:
                return self._button_factory.create(frame, **kwargs)
            elif widget_type == 'Label':
                kwargs.pop('type', None)
                return widget_set.Label(frame, **kwargs)
            else:
                raise NotImplementedError('type of widget not implemented')
