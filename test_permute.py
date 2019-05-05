import qutip as qt


base = qt.cnot()

state = qt.tensor(base, qt.qeye(2))

print(state)

plist = [0, 1, 2]
print(plist)
print(state.permute(plist))

plist = [1, 0, 2]
print(plist)
print(state.permute(plist))

plist = [0, 2, 1]
print(plist)
print(state.permute(plist))

plist = [1, 2, 0]
print(plist)
print(state.permute(plist))

plist = [2, 0, 1]
print(plist)
print(state.permute(plist))

plist = [2, 1, 0]
print(plist)
print(state.permute(plist))
