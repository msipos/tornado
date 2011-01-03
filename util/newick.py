import sys

class Consumer:
    def __init__(self, s):
        self.s = s
        self.i = 0
    
    def consume_char(self):
        c = self.s[self.i]
        self.i += 1
        return c

    def peek_char(self):
        c = self.s[self.i]
        return c

class Node:
    def __init__(self):
        self.children = []
    
    def write(self, f):
        l = len(self.children)
            
        if l > 0:
            f.write("(")
            i = 0
            for child in self.children:
                child.write(f, callback)
                if i < l - 1:
                    f.write(",")
                i += 1
            f.write(")")
        f.write(self.title)
        if "dist" in dir(self):
            f.write(":%f" % self.dist)

def write_newick(node, f):
    node.write(f)
    f.write(";")

def load_node(cons):
    n = Node()
    
    c = cons.peek_char()
    if c == "(":
        cons.consume_char()
        # Deal with children
        while True:
            node = load_node(cons)
            n.children.append(node)

            c = cons.consume_char()
            if c == ",":
                continue
            elif c == ")":
                break
            else:
                print "Error!"
    
    # Now, deal with title
    title = ""
    while True:
        c = cons.peek_char()
        if c != ":" and c!= ";":
            title += cons.consume_char()
        else:
            break
    n.title = title
    
    # Finally deal with distance
    c = cons.peek_char()
    if c == ":":
        cons.consume_char()
        text_dist = ""
        while True:
            c = cons.peek_char()
            if c in "0123456789.Ee-":
                text_dist += cons.consume_char()
            else:
                break
        dist = float(text_dist)
        n.dist = dist
    return n

def load_newick(filename):
    f = open(filename, "r")
    contents = f.read()
    cons = Consumer(contents)
    return load_node(cons)

