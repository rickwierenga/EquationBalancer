""" atoms.py """

import numpy as np
import sys
from numpy.linalg import det

# -- helpers
def gcd(x, y): 
    while(y): 
        x, y = y, x % y 

    return x


def gcd_list(a):
    result = a[0]
    for i in a[1:]:
        result = gcd(result, i)
    return result



# -- app
def molecules_in(string):
    """Get all the molecules in a string.
    Loop over the characters. If a character is a capital,
    add the current char to the dict. Else add the char to
    the current name.
    If the value is 0, 1 will be used as a value.
    if a value can be converted to an integer, add it to the
    current number.

    :returns: a dictionary of atoms with their number
    """
    strings = string.split(" + ")
    molecules = []

    for mol in strings:
        total_atoms = {}
        number_of_atoms = 0
        name = ''

        for char in mol:
            if char.isupper():
                if name != '':
                    total_atoms[name] = max(1, number_of_atoms)
                    number_of_atoms = 0
                    name = char
                else:
                    name += char
            else:
                try:
                    n = int(char)
                    number_of_atoms *= 10
                    number_of_atoms += n
                except ValueError:
                    name += char

        # add remainder
        total_atoms[name] = max(1, number_of_atoms)

        molecules.append(total_atoms)

    return molecules


def unique_atoms(molecules):
    """ Get all the unique atoms in an array of molecules
    """

    atoms = []
    for molecule in molecules:
        for atom in molecule.keys():
            if not atom in atoms:
                atoms.append(atom)
    return atoms


def inverse_values(molecules):
    inverse = []
    for mol in molecules:
        new_mol = {}
        for atom in mol.keys():
            new_mol[atom] = -1 * mol[atom]
        inverse.append(new_mol)
    return inverse


def construct_matrices(atoms, reactant, product):
    A_molecules = reactant + inverse_values(product[:-1]) 
    
    row_padding = 0
    column_padding = []
    if len(atoms) > len(A_molecules):
        row_padding = len(atoms) - len(A_molecules)
    elif len(A_molecules) > len(atoms):
        column_padding = [1] * (len(A_molecules) - len(atoms)) 

    A, B = [], []
    for atom in atoms:
        row = []
        for mol in A_molecules:
            try:
                row.append(mol[atom])
            except KeyError:
                row.append(0)
        row += [1] * row_padding
        A.append(row)
        
        try:
            B.append([product[-1][atom]])
        except KeyError:
            B.append([0])

    return np.array(A), np.array(B)


def normalise_answer(C):
    gcd_ = gcd_list(C)
    answer = []
    for num in C:
        answer.append(int(abs(num / gcd_)))
    return answer


def solve(reactant, product, atoms, A, B):
    A_ = np.linalg.inv(A)
    C = np.dot(A_, B) * det(A)
    print("it yields \n", C)

    # The number of molecules in A (= total minus the ones in B). The last item the determinant of A.
    num_mol_A = len(reactant) + len(product) - 1
    numbers = normalise_answer(list(C.reshape(1, len(C))[0][0:num_mol_A]) + [det(A)])

    reactant_answers = []
    product_answers = []

    for i in range(len(numbers)):
        try:
            reactant_answers.append([numbers[i], reactant[i]])
        except IndexError:
            product_answers.append([numbers[i], product[i - len(reactant)]])

    return reactant_answers, product_answers


def string_from_mol(idx_mol):
    idx = idx_mol[0]
    mol = idx_mol[1]

    if idx != 1:
        s = str(idx)
    else:
        s = ''

    for atom in mol.keys():
        if mol[atom] == 1:
            s += atom 
        else:
            s += atom + str(mol[atom])

    return s


def main():
    reaction = input("Please enter a reaction: ")
    try:
        # Get the molecules
        reactant_molecules, product_molecules = molecules_in(reaction.split(" -> ")[0]), molecules_in(reaction.split(" -> ")[1])
        print("I understood: ", reactant_molecules, " -> ", product_molecules)

        # Get the unique atoms
        atoms = unique_atoms(reactant_molecules)
        if not set(atoms) == set(unique_atoms(product_molecules)):
            print("reactant is not product")
            sys.exit(1)
        print("It contains the following atoms: ", atoms)
        
        # Create a one dimensional B matrix. Move the other molecules to A and multiply the number by -1.
        A, B = construct_matrices(atoms, reactant_molecules, product_molecules)
        print("Create 2 matrices: \n", A, "\n",  B)

        # Solve the matrix equation
        reactant_answer, product_answer = solve(reactant_molecules, product_molecules, atoms, A, B)

        # Create a string from the solution
        final = ''
        for mol in reactant_answer:
            final += string_from_mol(mol)
            final += " + "
        final = final[:-3] # remove superfluous +
        final += " -> "
        for mol in product_answer:
            final += string_from_mol(mol)
            final += " + "
        final = final[:-3] # remove superfluous +
        print(final)
    except IndexError:
        print("Invalid reaction")
        sys.exit(1)


if __name__ == '__main__':
    main()

