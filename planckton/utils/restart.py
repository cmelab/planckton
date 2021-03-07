import inspect

def _match_class_path(obj, *matches):
     return any(cls.__module__ + "." + cls.__name__ in matches
             for cls in inspect.getmro(type(obj)))
