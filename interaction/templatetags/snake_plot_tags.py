from django import template

register = template.Library()

@register.simple_tag
def render_snake_plot(structure, chain, nobuttons=None):
    print(structure, chain, nobuttons)
    return structure.get_snake_plot(chain=chain, nobuttons=nobuttons)
