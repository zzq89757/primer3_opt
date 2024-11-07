from collections import defaultdict, deque


int_loop_energy_dict = defaultdict(deque)
int_loop_energy_dict["22dh"] = deque([
		# A A
		# T T
		[
			6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
			7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
			7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
			6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
			7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
			7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
			6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
			7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
			7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
			6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
			7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
			7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
			6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
		],
		# A C
		# T G
		[
			-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
			-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
			-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
			-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
			-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
			-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
			-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
			-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
			-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
			-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
			-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
			-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
			-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
		],
		# A G
		# T C
		[
			-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
			0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
			0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
			-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
			0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
			0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
			-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
			0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
			0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
			-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
			0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
			0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
			-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
		],
		# A T
		# T A
		[
			8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
			9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
			8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
			8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
			8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
			9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
			8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
			8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
			8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
			9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
			8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
			8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
			8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
			9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
			8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
			8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
		],
		# C A
		# G T
		[
			-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
			-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
			-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
			-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
			-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
			-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
			-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
			-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
			-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
			-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
			-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
			-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
			-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
		],
		# C C
		# G G
		[
			-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
			-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
			-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
			-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
			-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
			-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
			-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
			-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
			-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
			-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
			-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
			-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
			-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
		],
		# C G
		# G C
		[
			-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
			-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
			-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
			-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
			-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
			-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
			-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
			-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
			-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
			-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
			-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
			-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
			-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
			-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
			-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
			-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
		],
		# C T
		# G A
		[
			-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
			-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
			-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
			-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
			-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
			-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
			-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
			-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
			-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
			-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
			-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
			-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
			-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
			-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
		],
		# G A
		# C T
		[
			-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
			-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
			-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
			1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
			-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
			-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
			1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
			-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
			-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
			1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
			-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
			-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
			1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
		],
		# G C
		# C G
		[
			-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
			-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
			-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
			-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
			-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
			-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
			-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
			-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
			-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
			-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
			-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
			-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
			-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
			-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
			-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
			-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
		],
		# G G
		# C C
		[
			-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
			-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
			-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
			-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
			-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
			-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
			-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
			-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
			-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
			-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
			-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
			-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
			-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
			-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
		],
		# G T
		# C A
		[
			-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
			0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
			-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
			3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
			0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
			-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
			3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
			0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
			-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
			3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
			-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
			0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
			-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
			3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
		],
		# T A
		# A T
		[
			4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
			4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
			3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
			-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
			4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
			4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
			3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
			-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
			4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
			4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
			3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
			-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
			4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
			4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
			3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
			-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
		],
		# T C
		# A G
		[
			-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
			-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
			-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
			-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
			-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
			-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
			-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
			-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
			-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
			-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
			-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
			-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
			-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
			-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
		],
		# T G
		# A C
		[
			-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
			-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
			-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
			-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
			-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
			-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
			-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
			-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
			-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
			-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
			-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
			-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
			-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
			-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
		],
		# T T
		# A A
		[
			6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
			6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
			4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
			-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
			6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
			4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
			-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
			6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
			4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
			-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
			6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
			6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
			4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
			-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
		],
])



matrix = int_loop_energy_dict["22dh"]
blank_li = []
an_averages = []
cn_averages = []
gn_averages = []
tn_averages = []
na_averages = []
nc_averages = []
ng_averages = []
nt_averages = []
nn_averages = []
num_columns = len(matrix[0])

for col in  range(num_columns):
    # AN (0-3)
    col_sum = sum(row[col] for i, row in enumerate(matrix) if i <= 3)  # 计算该列的总和
    an_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    an_averages.append(an_col_avg)  # 将平均值添加到结果列表
    # CN (4-7)
    col_sum = sum(row[col] for i, row in enumerate(matrix) if 3 < i <= 7)  # 计算该列的总和
    cn_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    cn_averages.append(cn_col_avg)  # 将平均值添加到结果列表
    # GN (8-11)
    col_sum = sum(row[col] for i, row in enumerate(matrix) if 7 < i <= 11)  # 计算该列的总和
    gn_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    gn_averages.append(gn_col_avg)  # 将平均值添加到结果列表
    # TN (12-15)
    col_sum = sum(row[col] for i, row in enumerate(matrix) if 11 < i <= 15)  # 计算该列的总和
    tn_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    tn_averages.append(tn_col_avg)  # 将平均值添加到结果列表
    
    
    # NA
    col_sum = sum(row[col] for i, row in enumerate(matrix) if i % 4 == 0)  # 计算该列的总和
    na_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    na_averages.append(na_col_avg)  # 将平均值添加到结果列表
    # NC
    col_sum = sum(row[col] for i, row in enumerate(matrix) if i % 4 == 1)  # 计算该列的总和
    nc_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    nc_averages.append(nc_col_avg)  # 将平均值添加到结果列表
    # NG
    col_sum = sum(row[col] for i, row in enumerate(matrix) if i % 4 == 2)  # 计算该列的总和
    ng_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    ng_averages.append(ng_col_avg)  # 将平均值添加到结果列表
    # NT
    col_sum = sum(row[col] for i, row in enumerate(matrix) if i % 4 == 3)  # 计算该列的总和
    nt_col_avg = round(col_sum / 4, 1)  # 计算该列的平均值
    nt_averages.append(nt_col_avg)  # 将平均值添加到结果列表
    
    # NN
    col_sum = sum(row[col] for row in matrix)  # 计算该列的总和
    nn_col_avg = round(col_sum / len(matrix), 1)  # 计算该列的平均值
    nn_averages.append(nn_col_avg)  # 将平均值添加到结果列表

# print(an_averages)

# insert (ACGT)N
inserted_cols = 0

int_loop_energy_dict["22dh"].insert(4,an_averages)
inserted_cols += 1

int_loop_energy_dict["22dh"].insert((4 + inserted_cols * 5),cn_averages)
inserted_cols += 1

int_loop_energy_dict["22dh"].insert((4 + inserted_cols * 5),gn_averages)
inserted_cols += 1

int_loop_energy_dict["22dh"].insert((4 + inserted_cols * 5),tn_averages)
inserted_cols += 1

# insert N(ACGTN)
int_loop_energy_dict["22dh"].append(na_averages)
int_loop_energy_dict["22dh"].append(nc_averages)
int_loop_energy_dict["22dh"].append(ng_averages)
int_loop_energy_dict["22dh"].append(nt_averages)
int_loop_energy_dict["22dh"].append(nn_averages)

# print filled matrix
basen_li = ["A", "C", "G", "T", "N"]

for mx_no, mx in enumerate(int_loop_energy_dict["22dh"]):
    first_base_idx = mx_no // 5
    second_base_idx = mx_no % 5
    print(f"\t\t# {basen_li[first_base_idx]}{basen_li[second_base_idx]}")
    print("\t\t[")
    # print(", ".join(mx))
    for idx, ele in enumerate(mx):
        if idx % 16 == 0:
            print("\t\t\t",end="")
        print(ele,end=", ")
        if idx % 16 == 15:
            print()
    print("\t\t],")

# print(len(int_loop_energy_dict["22dh"]))