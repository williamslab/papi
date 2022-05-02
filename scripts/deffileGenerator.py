def generateLines(gens):
    """
    Args:
    gens:Number of generations to simulate. 1 + observed number of crossovers
    """
#{{{ 
    lines = ['def {}gens 1250 {}'.format(gens-1,gens+1)]  # Header line. gens-1 is used so that header is consistent with observed crossover number
    remainingFounders = 2**(gens-1)
    current_gen = 1 

    if gens == 1:
        print('\n'.join(['def 0gens 1250 2','1 1 1','2 1 1'])) #Base case
        return()
    
    if gens == 2:
        print('\n'.join(['def 1gens 1250 3','1 0 2','2 1 2','3 1 1 1:1_2']))
        return()

    while remainingFounders > 1:

        if current_gen == 1 or current_gen == 2:
            current_line = ' '.join([str(x) for x in [current_gen,0,remainingFounders]])
        else:
            remainingFounders = remainingFounders//2

            num_print = 0 if remainingFounders > 2 else 1
            current_line_header = ' '.join([str(x) for x in [current_gen,num_print,remainingFounders]])
            current_line_branches = ' '.join(['{}:{}_{}'.format(i, i*2 - 1, i*2) for i in range(1,remainingFounders+1)])

            current_line = current_line_header + ' ' + current_line_branches

        lines.append(current_line)
        current_gen += 1

    print('\n'.join(lines))
#}}}
    return()

if __name__ == "__main__":
    import sys 
    generateLines(int(sys.argv[1]))
