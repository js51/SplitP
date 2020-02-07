def newick_to_json(newick_string, namestring = "id", lengthstring = "branch_length", childrenstring = "children", generate_names=False):
	"""Parses a newick string (example "((A,B),C);") into JSON format accepted by NetworkX

	Args:
		newick_string: The newick string to convert to JSON. Names, branch lengths and named 
			internal nodes accepted. If no names are given, an ID is generated for the name if generate_names=True.
		namestring: The label for node name to use in the JSON (default "id").
		lengthstring: The label for branch length to use in the JSON (default "branch_length").
		childrenstring: The label for children arrays to use in the JSON (default "children").
	"""
	def splitIntoChildren(children_string):
		""" Helper function for splitting a newick string into siblings """
		string = children_string
		opened = 0
		closed = 0
		pieces = []
		last_stop = 0
		for i, s in enumerate(string):
			if s == "(":
				opened += 1
			if s == ")":
				closed += 1
			if s == "," and opened == closed:
				pieces.append(string[last_stop:i])
				last_stop = i+1
		pieces.append(string[last_stop:])
		return(pieces)
		
	try: 
		newick_to_json.counter += 1
	except:
		newick_to_json.counter = 0
	import re
	newick_string = newick_string.strip(";")
	node = dict()
	info = newick_string[newick_string.rfind(")")+1:]
	info = info.split(":")
	name = info[0]
	if name == "" and generate_names:
		node[namestring] = newick_to_json.counter
	else:
		node[namestring] = name
	length = ""
	if len(info) > 1:
		length = info[1]
		node[lengthstring] = float(length)
	children_string = newick_string[0:newick_string.rfind(")")+1]
	if children_string != "":
		if children_string[-1] == ")":
				children_string = children_string[1:-1]
		children_string = splitIntoChildren(children_string)
		child_nodes_JSON = []
		for child_string in children_string:
			child_nodes_JSON.append(newick_to_json(child_string, namestring, lengthstring, childrenstring))
		node[childrenstring] = child_nodes_JSON
	return node
	