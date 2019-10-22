"""Simple analysis of a cantilevered pipe"""

# inputs
L = 10  # ft

Model('simple')

# attributes
Pipe.from_file('pipe1', '10', '40')
Material.from_file('mat1', 'A53A', 'B31.1')

# geometry
Point(10)
Run(20, L)

# supports
Anchor('A1', 10)

# loads
W1 = Weight('W1', 386.6)
P1 = Pressure('P1', 100)
elements.apply_loads(W1, P1)

# loadcases
L1 = LoadCase('L1', [W1, P1], "OPE")
L1.solve()

# results
disp = results.Displacement(L1)
disp.to_text('disp.out')

app.quit()
