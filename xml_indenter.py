
def smart_indent_xml(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            smart_indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
    return elem

def split_attributes(line):
    start = line.index('<')
    indent = line[:start]
    end = line.index('>')
    between = line[start+1:end]
    parts = between.split(' ')
    tag = parts[0]
    attributes = {}
    for part in parts[1:]:
        key, value = part.split('=')
        attributes[key] = value
    return indent, tag, attributes

def new_split_line(tag, attributes, indent, gap):
    line = indent + '<' + tag
    linenum = 0
    for key, value in attributes.items():
        extra = ' ' * (((gap+1)*(linenum>0)) + 1)
        line += indent + extra  + key + '=' + value + '\n'
        linenum += 1
    line = line.rstrip()  
    line += '>\n'
    return line

def split_xml_attributes(xmlfile, field = '<'):
    with open(xmlfile, 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if (line.lstrip()[0] == '<') and (field in line):
            indent, tag, attributes = split_attributes(line)
            new_line = new_split_line(tag, attributes, indent, len(tag))
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    with open(xmlfile, 'w') as f:
        f.writelines(new_lines)
    return new_lines
 