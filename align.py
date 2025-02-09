"""
Created on 1/29/25
@author: Mohammad Abbas Naqvigit
"""


#### linear pen == 1

"""
myscoring = {'AA' : 5, 'AG' :-1, 'GA': -1, 'AT':-4, 'TA':-4, 'AC':-4,
             'CA':-4, 'GG':5,'GC':-4, 'CG':-4, 'GT':-4,'TG':-4,'CC':5,
             'TC':-1,'CT':-1,'TT':5}
"""

myscoring = {'AA' : 5, 'AG' :-1, 'GA': -1, 'AT':-1, 'TA':-1, 'AC':-1,
             'CA':-1, 'GG':5,'GC':-1, 'CG':-1, 'GT':-1,'TG':-1,'CC':5,
             'TC':-1,'CT':-1,'TT':5}

def linear_gap(x, y):
    matrix = [[0 for each in range(len(x) + 1)] for each in range(len(y) + 1)]
    trace_matrix = [['' for each in range(len(x) + 1)] for each in range(len(y) + 1)]

    x_list = list(x)
    y_list = list(y)


    print(x_list, y_list)

### Initialize borders ###
    i = 0
    for every in range(len(x) + 1):
        matrix[0][every] += i
        i += -1

    i = 0
    for every in range(len(y) + 1):
        matrix[every][0] += i
        i += -1

### Fill in matrix ###
    for row in range(1, len(y) + 1):
        for column in range(1, len(x) + 1):
            diag_cell = matrix[row - 1][column - 1] + myscoring[x_list[column-1]+y_list[row-1]]
            left_cell = matrix[row][column-1] - 1
            up_cell = matrix[row-1][column] - 1

            matrix[row][column] = max(left_cell, up_cell,diag_cell)

            if(max(left_cell, up_cell,diag_cell)) == diag_cell:
                trace_matrix[row][column] = "Diagonal"
            if(max(left_cell, up_cell,diag_cell)) == left_cell:
                trace_matrix[row][column] = "Left"
            if(max(left_cell, up_cell,diag_cell)) == up_cell:
                trace_matrix[row][column] = "Up"





    for line in trace_matrix:
        print(line)
    for line in matrix:
        print(line)



#Optimal Alignment Code

    opt_X = []
    opt_Y = []
    seqX_length = len(x)
    seqY_length =  len(y)

    alignment_score = 0
    while seqX_length > 0 or seqY_length > 0:
        if trace_matrix[seqY_length][seqX_length] == 'Diagonal':
            opt_X.append(x_list[seqX_length-1])
            opt_Y.append(y_list[seqY_length-1])
            seqX_length -= 1
            seqY_length -= 1

        if trace_matrix[seqY_length][seqX_length] == 'Left':
            opt_X.append(x_list[seqX_length-1])
            opt_Y.append('-')
            seqX_length -= 1
        
        if trace_matrix[seqY_length][seqX_length] == 'Up':
            opt_X.append('-')
            opt_Y.append(y_list[seqY_length-1])
            seqY_length -= 1

    #To reverse order because of traceback  
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
    print("Percent identity: " + str((num / (len(x_list) + len(y_list))*100)) + '%')
    print("Indels: number=" + str(indelquant)+ ", mean length=" + str(mean_length))
    print("Alignment length: " + str(len(final_x)))


    
    for x in range(0, len(final_x), 60):
        if len(final_x) > x + 60:
            print(final_x[x:x+60])
            print(dashstar[x:x+60])
            print(final_y[x:x+60])
        else:
            print(final_x[x:])
            print(dashstar[x:])
            print(final_y[x:])
    #print(alignment_score)


   


        
        


        
        



    return






def affine_gap(x, y):
    matrix = [[0 for each in range(len(x) + 1)] for each in range(len(y) + 1)]
    trace_matrix = [['' for each in range(len(x) + 1)] for each in range(len(y) + 1)]

    x_list = list(x)
    y_list = list(y)

    print(x_list, y_list)

    # Initialize borders
    i = 0
    for every in range(len(x) + 1):
        matrix[0][every] += i
        i += -1

    i = 0
    for every in range(len(y) + 1):
        matrix[every][0] += i
        i += -1

    # Fill in matrix
    isDiag = True
    for row in range(1, len(y) + 1):
        for column in range(1, len(x) + 1):
            diag_cell = matrix[row - 1][column - 1] + myscoring[x_list[column - 1] + y_list[row - 1]]
            if isDiag:
                left_cell = matrix[row][column - 1] - 4
                up_cell = matrix[row - 1][column] - 4
            else:
                left_cell = matrix[row][column - 1] - 1
                up_cell = matrix[row - 1][column] - 1

            matrix[row][column] = max(left_cell, up_cell, diag_cell)
            isDiag = matrix[row][column] == diag_cell

            if max(left_cell, up_cell, diag_cell) == diag_cell:
                trace_matrix[row][column] = "Diagonal"
            if max(left_cell, up_cell, diag_cell) == left_cell:
                trace_matrix[row][column] = "Left"
            if max(left_cell, up_cell, diag_cell) == up_cell:
                trace_matrix[row][column] = "Up"

    for line in trace_matrix:
        print(line)
    for line in matrix:
        print(line)

    # Optimal Alignment Code
    opt_X = []
    opt_Y = []
    seqX_length = len(x)
    seqY_length = len(y)

    alignment_score = 0
    while seqX_length > 0 or seqY_length > 0:
        if trace_matrix[seqY_length][seqX_length] == 'Diagonal':
            opt_X.append(x_list[seqX_length - 1])
            opt_Y.append(y_list[seqY_length - 1])
            seqX_length -= 1
            seqY_length -= 1

        elif trace_matrix[seqY_length][seqX_length] == 'Left':
            opt_X.append(x_list[seqX_length - 1])
            opt_Y.append('-')
            seqX_length -= 1

        elif trace_matrix[seqY_length][seqX_length] == 'Up':
            opt_X.append('-')
            opt_Y.append(y_list[seqY_length - 1])
            seqY_length -= 1

    # Reverse order due to traceback
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
    print("Percent identity:", str((num / (len(x_list) + len(y_list)) * 100)) + '%')
    print("Indels: number =", str(indelquant) + ", mean length =", str(mean_length))
    print("Alignment length:", len(final_x))

    for x in range(0, len(final_x), 60):
        if len(final_x) > x + 60:
            print(final_x[x:x+60])
            print(dashstar[x:x+60])
            print(final_y[x:x+60])
        else:
            print(final_x[x:])
            print(dashstar[x:])
            print(final_y[x:])

    #print(alignment_score)




if __name__ == '__main__':
    print('hello')
    print(linear_gap('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA','GATAACATTGATACGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA'))
    print(affine_gap('GATTACA','GATAACA'))
