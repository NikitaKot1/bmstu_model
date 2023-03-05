import analit
import algos
from functio import *
from prettytable import PrettyTable



def task_1() -> PrettyTable:
    tabl = PrettyTable()

    tabl.field_names = ["Аргумент", "Аналит.", "Эйлер", "Пикар 1", "Пикар 2", "Пикар 3", "Пикар 4"]

    x_min = 1
    x_max = 2
    y_start = 0
    h = 1e-3
    n = int((x_max - x_min) / h)
    eul = algos.euler(x_min, y_start, h, n, f_1)
    picar1 = algos.picar(x_max, h, algos.picar_apox_11, 1, 0)
    picar2 = algos.picar(x_max, h, algos.picar_apox_12, 1, 0)
    picar3 = algos.picar(x_max, h, algos.picar_apox_13, 1, 0)
    picar4 = algos.picar(x_max, h, algos.picar_apox_14, 1, 0)
    k = -1
    for i in range(n):
        row = []
        row.append("{:^9.3f}".format(x_min))
        row.append("{:^9.5f}".format(analit.analit_1(x_min)))
        if eul[i] > 10000000:
            row.append("Too big")
        else:
            row.append("{:^9.5f}".format(eul[i]))

        row.append("{:^9.5f}".format(picar1[i]))
        row.append("{:^9.5f}".format(picar2[i]))
        row.append("{:^9.5f}".format(picar3[i]))
        row.append("{:^9.5f}".format(picar4[i]))
        k += 1
        x_min += h

        if k % 5 == 0:
            tabl.add_row(row)

    return tabl

def task_2() -> PrettyTable:
    tabl = PrettyTable()

    tabl.field_names = ["Аргумент", "Аналит.", "Эйлер", "Пикар 1", "Пикар 2", "Пикар 3", "Пикар 4"]

    x_min = 0.5
    x_max = 2
    y_start = 0
    h = 1e-2
    n = int((x_max - x_min) / h)
    eul = algos.euler(x_min, y_start, h, n, f_2)
    picar1 = algos.picar(x_max, h, algos.picar_apox_21, 0.5, 0)
    picar2 = algos.picar(x_max, h, algos.picar_apox_22, 0.5, 0)
    picar3 = algos.picar(x_max, h, algos.picar_apox_23, 0.5, 0)
    picar4 = algos.picar(x_max, h, algos.picar_apox_24, 0.5, 0)
    k = -1
    for i in range(n):
        row = []
        row.append("{:^9.3f}".format(x_min))
        row.append("{:^9.5f}".format(analit.analit_2(x_min)))
        if eul[i] > 10000000:
            row.append("Too big")
        else:
            row.append("{:^9.5f}".format(eul[i]))

        row.append("{:^9.5f}".format(picar1[i]))
        row.append("{:^9.5f}".format(picar2[i]))
        row.append("{:^9.5f}".format(picar3[i]))
        row.append("{:^9.5f}".format(picar4[i]))
        k += 1
        x_min += h

        if k % 10 == 0:
            tabl.add_row(row)


    return tabl

def task_3() -> PrettyTable:
    tabl = PrettyTable()

    tabl.field_names = ["Аргумент", "Эйлер", "Пикар 1", "Пикар 2", "Пикар 3", "Пикар 4"]

    x_min = 0
    x_max = 1.3
    y_start = 0
    h = 1e-3
    n = int((x_max - x_min) / h)
    eul = algos.euler(x_min, y_start, h, n, f_3)
    picar1 = algos.picar(x_max, h, algos.picar_apox_31, 0, 0)
    picar2 = algos.picar(x_max, h, algos.picar_apox_32, 0, 0)
    picar3 = algos.picar(x_max, h, algos.picar_apox_33, 0, 0)
    picar4 = algos.picar(x_max, h, algos.picar_apox_34, 0, 0)
    k = -1
    for i in range(n):
        row = []
        row.append("{:^9.3f}".format(x_min))
        if type(eul[i]) == int:
            if eul[i] > 10000000:
                row.append("Too big")
            else:
                row.append("{:^9.5f}".format(eul[i]))
        else:
            row.append(eul[i])

        row.append("{:^9.5f}".format(picar1[i]))
        row.append("{:^9.5f}".format(picar2[i]))
        row.append("{:^9.5f}".format(picar3[i]))
        row.append("{:^9.5f}".format(picar4[i]))
        k += 1
        x_min += h

        if k % 10 == 0:
            tabl.add_row(row)


    return tabl


def main():
    task_num = int(input("Введите номер задания: "))
    if task_num == 1:
        tabl = task_1()
    elif task_num == 2:
        tabl = task_2()
    else:
        tabl = task_3()
    print(tabl)


if __name__ == "__main__":
    main()



# for i in range(10):
#     print(i / 10, analit.analit(i / 10))