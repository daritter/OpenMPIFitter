import os
from lxml import etree

from SCons.Action import Action
from SCons.Builder import Builder

def scan_svg(target,source,env):
    svgfile = open(source[0].rfile().abspath,'r')
    svg = etree.parse(svgfile)
    for img in svg.xpath('//svg:image',namespaces={"svg":"http://www.w3.org/2000/svg"}):
        href =  img.get('{http://www.w3.org/1999/xlink}href')
        if href.startswith('data:'): continue
        if href.startswith('file://'): href = href[7:]
        if not os.path.isabs(href):
            href = os.path.join(os.path.dirname(source[0].rfile().abspath),href)

        source.append(os.path.normpath(href))

    target.append(str(target[0])+"_tex")

    return target,source

build_svg = "inkscape -z -C --export-latex --export-pdf=$TARGET $SOURCE"

buildSvgAction = Action(build_svg,'${_printCMD("purple","ink2pdftex")} $TARGET')

ink2pdf = Builder(
    action=buildSvgAction,
    suffix='.pdf',
    src_suffix='.svg',
    emitter=scan_svg,
)

def generate(env):
    env.Tool('output')
    env['BUILDERS']['InkPDF'] = ink2pdf

def exists(env):
    return True
