from django import template

register = template.Library()

@register.filter(name="get_index")
def get_index(sequence, index):
    """Return item at numeric index from a list/tuple; None if out of range."""
    try:
        return sequence[index]
    except Exception:
        return None
