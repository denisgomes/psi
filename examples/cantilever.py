"""Simple analysis of a cantilevered pipe"""

# inputs
L = 10

mdl = Model('simple')

# attributes
pipe1 = Pipe.from_file('pipe1', '10', '40')
mat1 = Material.from_file('mat1', 'A53A', 'B31.1')

# geometry
pt10 = Point(10)
run20 = Run(20, L)

# supports
Anchor('A1', 10)

# loads
f1 = Force('f1', 20, fy=1000)
f1.apply()  # to active element

# loadcases
l1 = LoadCase('L1', [f1], 'OPE')
l1.solve()

# results
# disp = results.Displacement(l1)
# disp.to_text('disp.out')

# app.quit()
