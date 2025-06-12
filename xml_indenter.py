import shlex


# print out root and children as a tree showing attributes and values, recursively
def print_tree(node, level=0, print_attributes=False):
    # print('-'*20)
    print("  " * level + node.tag)
    # print the attributes
    if print_attributes:
        for key, value in node.attrib.items():
            print("  " * (level + 1) + key, ":", value)
    for child in node:
        print_tree(child, level + 1, print_attributes=print_attributes)


def smart_indent_xml(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            smart_indent_xml(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
    return elem


def split_attributes(line):
    start = line.index("<")
    indent = line[:start]
    end = line.index(">")
    between = line[start + 1 : end]
    parts = shlex.split(between)
    tag = parts[0]
    attributes = {}
    for part in parts[1:]:
        key, value = part.split("=")
        value = '"' + value + '"'
        value.replace('""', '"')
        attributes[key] = value
    return indent, tag, attributes


def new_split_line(tag, attributes, indent, gap):
    line = indent + "<" + tag
    linenum = 0
    for key, value in attributes.items():
        extra = " " * (((gap + 1) * (linenum > 0)) + 1)
        line += indent + extra + key + "=" + value + "\n"
        linenum += 1
    line = line.rstrip()
    line += ">\n"
    return line


def split_xml_attributes(xmlfile, field="<"):
    with open(xmlfile, "r") as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if (line.lstrip()[0] == "<") and (field in line):
            indent, tag, attributes = split_attributes(line)
            new_line = new_split_line(tag, attributes, indent, len(tag))
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    with open(xmlfile, "w") as f:
        f.writelines(new_lines)
    return new_lines


def replace_lt_gt(xmlfile):
    with open(xmlfile, "r") as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if line.lstrip()[0] == "<":
            line = line.replace("&lt;", "<")
            line = line.replace("&gt;", ">")
        new_lines.append(line)
    with open(xmlfile, "w") as f:
        f.writelines(new_lines)
    return new_lines
