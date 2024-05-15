n = 5
mod = 31
# generator for my lib

print("std::vector<Monomial<PrimeField<" + str(mod) + ">>> monk;")
print()

print("monk.push_back(Monomial(Term({1}), PrimeField<" + str(mod) + ">(1)));")

tt = "0, "
for i in range(1, n):
    kk = tt + "1"
    tt += "0, "
    print("monk.push_back(Monomial(Term({" + kk + "}), PrimeField<" + str(mod) + ">(2)));")
print("monk.push_back(Monomial(Term({0}), PrimeField<" + str(mod) + ">(-1)));")
print()
print("Polynomial<PrimeField<" + str(mod) + ">, GrevLexComp> pk(std::move(monk));")
print()

for i in range(n - 1):
    print("std::vector<Monomial<PrimeField<" + str(mod) + ">>> mon" + str(i) + ";")
    print()
    arr = []
    for j in range(n):
        if abs(i - j) > j:
            continue

        ar = [0] * (j + 1)
        cf = 2
        if i / 2 != j:
            ar[j] += 1
            ar[(abs(i - j))] += 1
        else:
            ar = [0] * (j + 1)
            ar[j] = 2
            cf = 1
        arr.append([ar, cf])
    ar = [0] * (i + 1)
    ar[i] = 1
    arr.append([ar, -1])
    for x in arr:
        t = ""
        for k in range(len(x[0])):
            t += str(x[0][k])
            if k != len(x[0]) - 1:
                t += ", "
        print("mon" + str(i) + ".push_back(Monomial(Term({" + t + "}), PrimeField<" + str(mod) + ">(" + str(x[1]) + ")));")
    print()
    print("Polynomial<PrimeField<" + str(mod) + ">, GrevLexComp> p" + str(i) + "(std::move(mon" + str(i) + "));")
    print()

finish = "pk, "
for i in range(n - 1):
    finish += "p" + str(i)
    if i != n - 2:
        finish += ", "
print("TPolynomials<PrimeField<" + str(mod) + ">, GrevLexComp> test = {" + finish + "};")
