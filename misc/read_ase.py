from ase.db import connect

with connect('Si_ase.db') as db:
    for i, row in enumerate(db.select()):
        print(row.data['energy'])

