"""
Created on 1/29/25
@author: Mohammad Abbas Naqvi
"""


#### linear pen == 1

def parser(f_name):
	sequence_Dict = {}
	with open(f_name, 'r') as file:
		curr_seq = []
		curr_tag = 0
		for line in file:
			myline = line.strip()
			if myline and myline[0] == ">":
				if curr_tag != 0:
					sequence_Dict[curr_tag] = ''.join(curr_seq)
				curr_tag = myline[1:]
				curr_seq = []
			elif myline:
				curr_seq.append(myline)
		if curr_tag != 0:
			sequence_Dict[curr_tag] = ''.join(curr_seq)
	return sequence_Dict


myscoring = {'AA' : 5, 'AG' :-1, 'GA': -1, 'AT':-4, 'TA':-4, 'AC':-4,
             'CA':-4, 'GG':5,'GC':-4, 'CG':-4, 'GT':-4,'TG':-4,'CC':5,
             'TC':-1,'CT':-1,'TT':5}
			 

"""
myscoring = {'AA' : 5, 'AG' :-1, 'GA': -1, 'AT':-1, 'TA':-1, 'AC':-1,
             'CA':-1, 'GG':5,'GC':-1, 'CG':-1, 'GT':-1,'TG':-1,'CC':5,
             'TC':-1,'CT':-1,'TT':5}
			 """
             

def linear_gap(x, y):
    matrix = [[0 for each in range(len(x) + 1)] for each in range(len(y) + 1)]
    trace_matrix = [['' for each in range(len(x) + 1)] for each in range(len(y) + 1)]

    x_list = list(x)
    y_list = list(y)



### Initialize borders ###
    
    i = 0
    for every in range(len(x) + 1):
        matrix[0][every] += i
        i += -1

    i = 0
    for every in range(len(y) + 1):
        matrix[every][0] += i
        i += -1

    for i in range(1, len(x) + 1):
        trace_matrix[0][i] = 'Left'


    for i in range(1, len(y) + 1):
        trace_matrix[i][0] = 'Up'

		

### Fill in matrix ###
    for row in range(1, len(y) + 1):
        for column in range(1, len(x) + 1):
            diag_cell = matrix[row - 1][column - 1] + myscoring[x_list[column-1]+y_list[row-1]]
            left_cell = matrix[row][column-1] - 1
            up_cell = matrix[row-1][column] - 1

            matrix[row][column] = max(left_cell, up_cell,diag_cell)

            if(max(left_cell, up_cell,diag_cell)) == diag_cell:
                trace_matrix[row][column] = "Diagonal"
            elif(max(left_cell, up_cell,diag_cell)) == left_cell:
                trace_matrix[row][column] = "Left"
            else:
                trace_matrix[row][column] = "Up"


    

#Optimal Alignment Code

    opt_X = []
    opt_Y = []
    seqX_length = len(x)
    seqY_length =  len(y)

    alignment_score = matrix[len(y)][len(x)]

    
    while seqX_length > 0 or seqY_length > 0:
        if trace_matrix[seqY_length][seqX_length] == 'Diagonal':
            opt_X.append(x_list[seqX_length-1])
            opt_Y.append(y_list[seqY_length-1])
            seqX_length -= 1
            seqY_length -= 1

        elif trace_matrix[seqY_length][seqX_length] == 'Left':
            opt_X.append(x_list[seqX_length-1])
            opt_Y.append('-')
            seqX_length -= 1
        
        elif trace_matrix[seqY_length][seqX_length] == 'Up':
            opt_X.append('-')
            opt_Y.append(y_list[seqY_length-1])
            seqY_length -= 1
        else:
            print("whoops")
        
   
        
    opt_X.reverse()
    opt_Y.reverse()

    dashstar_list = ['' for each in range(len(opt_X))]

    for i in range(0,len(opt_X)):
        if opt_X[i] == opt_Y[i]:
            dashstar_list[i] = '|'
        elif (opt_X[i] != opt_Y[i]) and ((opt_X[i] == '-') or (opt_Y[i] == '-')):
            dashstar_list[i] = ' '
        else:
            dashstar_list[i] = '*'
        
    final_x = ''.join(opt_X)
    dashstar = ''.join(dashstar_list)
    final_y = ''.join(opt_Y)


    indel = []
    indelquant = 0
  

    i = 0
    while i < len(opt_X):
        if opt_X[i] == '-':
            count = 0
            while i < len(opt_X) and opt_X[i] == '-':
                count += 1
                i += 1
            indel.append(count)  
            indelquant += 1  
        else:
            i += 1

    i = 0
    while i < len(opt_Y):
        if opt_Y[i] == '-':
            count = 0
            while i < len(opt_Y) and opt_Y[i] == '-':
                count += 1
                i += 1
            indel.append(count)  
            indelquant += 1  
        else:
            i += 1

    length = 0
    for i in indel:
        length += i
    
    if len(indel) != 0:
        mean_length = str(length / len(indel))
    else:
        mean_length = 0
    
    


    num = dashstar.count('|')
    print("Matches: " + str(num))
    print("Percent identity: " + str((num / ((len(x_list) + len(y_list)) / 2)) * 100) + '%')
    print("Indels: number=" + str(indelquant)+ ", mean length=" + str(mean_length))
    print("Alignment length: " + str(len(final_x)))
    print("Score=" + str(alignment_score))



    
    for x in range(0, len(final_x), 60):
        if len(final_x) > x + 60:
            print(final_x[x:x+60])
            print(dashstar[x:x+60])
            print(final_y[x:x+60])
        else:
            print(final_x[x:])
            print(dashstar[x:])
            print(final_y[x:])



    return






def affine_gap(x, y):

    Match_matrix = [[float(0) for each in range(len(x) + 1)] for each in range(
        len(y) + 1)]
    Insertion_matrix = [[float(0) for each in range(len(x) + 1)] for each in range(
        len(y) + 1)]

    Deletion_matrix = [[float(0) for each in range(len(x) + 1)] for each in
                    range(len(y) + 1)]

    T_match_matrix = [['' for each in range(len(x) + 1)] for each in
                    range(len(y) + 1)]

    T_insertion_matrix = [['' for each in range(len(x) + 1)] for each in
                    range(len(y) + 1)]

    T_deletion_matrix = [['' for each in range(len(x) + 1)] for each in
                    range(len(y) + 1)]

    x_list = list(x)
    y_list = list(y)

    m = len(x_list)
    n = len(y_list)

    Deletion_matrix[0][0] = low
    Insertion_matrix[0][0] = low

    for i in range(1, m + 1):
        Match_matrix[0][i] = low
        Deletion_matrix[0][i] = low


        Insertion_matrix[0][i] = max(Insertion_matrix[0][i-1] - 0.1,
                                     Match_matrix[0][i-1] - 4)
        if Insertion_matrix[0][i] == Insertion_matrix[0][i-1] - 0.1:
            T_insertion_matrix[0][i] = "Insertion"
        else:
            T_insertion_matrix[0][i] = "Match"

    for i in range(1, n + 1):
        Match_matrix[i][0] = low
        Insertion_matrix[i][0] = low

        Deletion_matrix[i][0] = max(Deletion_matrix[i-1][0] - 0.1,
                                     Match_matrix[i-1][0] - 4)
        if Deletion_matrix[i][0] == Deletion_matrix[i-1][0] - 0.1:
            T_deletion_matrix[i][0] = "Deletion"
        else:
            T_deletion_matrix[i][0] = "Match"

    for row in range(1,n+1):
        for col in range(1,m+1):
            
            #Recurrence for Match code below
            Match = Match_matrix[row-1][col-1] + myscoring[x_list[col-1] + y_list[row-1]]
            Insertion = Insertion_matrix[row-1][col-1] + myscoring[x_list[col-1] + y_list[row-1]]
            Deletion = Deletion_matrix[row-1][col-1] + myscoring[x_list[col-1] + y_list[row-1]]
            Match_matrix[row][col] = max(Match,Insertion,Deletion)
            if Match_matrix[row][col] == Match:
                T_match_matrix[row][col] = 'Match'
            elif Match_matrix[row][col] == Deletion:
                T_match_matrix[row][col] = "Deletion"
            else:
                T_match_matrix[row][col] = "Insertion"

            # Recurrence for Deletion Code Below
            extension_cum = Deletion_matrix[row-1][col] - 0.1
            open_cum = Match_matrix[row-1][col] - 4
            Deletion_matrix[row][col] = max(extension_cum, open_cum)
            if Deletion_matrix[row][col] == open_cum:
                T_deletion_matrix[row][col] = "Match"
            else:
                T_deletion_matrix[row][col] = "Deletion"

            #Recurrence for Insertion Code Below
            extension_cum = Insertion_matrix[row][col - 1] - 0.1
            open_cum = Match_matrix[row][col - 1] - 4
            Insertion_matrix[row][col] = max(extension_cum,open_cum)
            if Insertion_matrix[row][col] == open_cum:
                T_insertion_matrix[row][col] = "Match"
            else:
                T_insertion_matrix[row][col] = "Insertion"

    Match_score = Match_matrix[n][m]
    Insertion_score = Insertion_matrix[n][m]
    Deletion_score = Deletion_matrix[n][m]
    best_matrix = 'Match'
    best_score = Match_score
    if Insertion_score > best_score:
        best_matrix = 'Insertion'
        best_score = Insertion_score
    if Deletion_score > best_score:
        best_matrix = 'Deletion'
        best_score = Deletion_score

#m = len(x) = j = seqX
#n = len(y) = i = seqY

    seqX_length = m
    seqY_length = n
    matrix_temp = best_matrix
    opt_X = []
    opt_Y = []
    while seqX_length > 0 or seqY_length > 0:
        if matrix_temp == "Match":
            matrix_temp = T_match_matrix[seqY_length][seqX_length]
            opt_X.append(x_list[seqX_length - 1])
            opt_Y.append(y_list[seqY_length - 1])
            seqX_length -= 1
            seqY_length -= 1
        elif matrix_temp == 'Insertion':
            matrix_temp = T_insertion_matrix[seqY_length][seqX_length]
            opt_X.append(x_list[seqX_length - 1])
            opt_Y.append('-')
            seqX_length -= 1
        elif matrix_temp == 'Deletion':
            matrix_temp = T_deletion_matrix[seqY_length][seqX_length]
            opt_X.append('-')
            opt_Y.append(y_list[seqY_length - 1])
            seqY_length -= 1
        else:
            print("whoops")


    opt_X.reverse()
    opt_Y.reverse()

    dashstar_list = ['' for _ in range(len(opt_X))]

    for i in range(len(opt_X)):
        if opt_X[i] == opt_Y[i]:
            dashstar_list[i] = '|'
        elif opt_X[i] != opt_Y[i] and ('-' in (opt_X[i], opt_Y[i])):
            dashstar_list[i] = ' '
        else:
            dashstar_list[i] = '*'

    final_x = ''.join(opt_X)
    dashstar = ''.join(dashstar_list)
    final_y = ''.join(opt_Y)

    indel = []
    indelquant = 0

    i = 0
    while i < len(opt_X):
        if opt_X[i] == '-':
            count = 0
            while i < len(opt_X) and opt_X[i] == '-':
                count += 1
                i += 1
            indel.append(count)
            indelquant += 1
        else:
            i += 1

    i = 0
    while i < len(opt_Y):
        if opt_Y[i] == '-':
            count = 0
            while i < len(opt_Y) and opt_Y[i] == '-':
                count += 1
                i += 1
            indel.append(count)
            indelquant += 1
        else:
            i += 1

    length = sum(indel)
    if len(indel) != 0:
        mean_length = str(length / len(indel))
    else:
        mean_length = 0

    num = dashstar.count('|')
    print("Matches:", num)
    print("Percent identity:", str((num / ((len(x_list) + len(y_list)) / 2)) * 100) + '%')
    print("Indels: number =", str(indelquant) + ", mean length =", str(mean_length))
    print("Alignment length:", len(final_x))
    print("Score=" + str(best_score))
    print("\n")

    for x in range(0, len(final_x), 60):
        if len(final_x) > x + 60:
            print(final_x[x:x+60])
            print(dashstar[x:x+60])
            print(final_y[x:x+60])
        else:
            print(final_x[x:])
            print(dashstar[x:])
            print(final_y[x:])


if __name__ == '__main__':

    
    close_first_dict = parser("/Users/abbas/Desktop/Alignment/close-first.fasta" )
    close_second_dict = parser("/Users/abbas/Desktop/Alignment/close-second.fasta")
    distant_first_dict = parser("/Users/abbas/Desktop/Alignment/distant-first.fasta")
    distant_second_dict = parser("/Users/abbas/Desktop/Alignment/distant-second.fasta")

    close_first_key_list = [x for x in close_first_dict.keys()]
    close_second_key_list = [x for x in close_second_dict.keys()]
    distant_first_key_list = [x for x in distant_first_dict.keys()]
    distant_second_key_list = [x for x in distant_second_dict.keys()]

    print("\n")
    print("Distant Linear:\n")
    for i in range(0, len(close_first_key_list)):
        print("Alignment #" + str(i + 1) + ":\n\n")
        print("Sequence #1: " + distant_first_key_list[i])
        print("Sequence #2: " + distant_second_key_list[i])
        linear_gap(distant_first_dict[distant_first_key_list[i]], distant_second_dict[distant_second_key_list[i]])
        print("\n")


"""
    print("\n")
    print("Close Linear:\n")
    for i in range(0, len(close_first_key_list)):
        print("Alignment #" + str(i + 1) + ":\n\n")
        print("Sequence #1: " + close_first_key_list[i])
        print("Sequence #2: " + close_second_key_list[i])
        linear_gap(close_first_dict[close_first_key_list[i]], close_second_dict[close_second_key_list[i]])
        print("\n")
    

    print("\n")
    print("Close Affine:\n")
    for i in range(0, len(close_first_key_list)):
        print("Alignment #" + str(i + 1) + ":\n\n")
        print("Sequence #1: " + close_first_key_list[i])
        print("Sequence #2: " + close_second_key_list[i])
        affine_gap(close_first_dict[close_first_key_list[i]], close_second_dict[close_second_key_list[i]])
        print("\n")
        
    
    print("\n")
    print("Distant Affine:\n")
    for i in range(0, len(close_first_key_list)):
        print("Alignment #" + str(i + 1) + ":\n\n")
        print("Sequence #1: " + distant_first_key_list[i])
        print("Sequence #2: " + distant_second_key_list[i])
        affine_gap(distant_first_dict[distant_first_key_list[i]], distant_second_dict[distant_second_key_list[i]])
        print("\n")

        
    
    print("\n")
    print("Distant Linear:\n")
    for i in range(0, len(close_first_key_list)):
        print("Alignment #" + str(i + 1) + ":\n\n")
        print("Sequence #1: " + distant_first_key_list[i])
        print("Sequence #2: " + distant_second_key_list[i])
        linear_gap(distant_first_dict[distant_first_key_list[i]], distant_second_dict[distant_second_key_list[i]])
        print("\n")
        """

        
        
        
    
        

     

    


 


    

    

    


    #print('hello')
    #print(linear_gap('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA','GATAACATTGATACGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA'))
    #print(affine_gap('GACGAGAGTAGCGTATGAGCGTA','GATGATCGTAGTCGAGTGCTGAT'))
    


