# Добавление зависимостей для работы с FastAPI
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
# Добавление зависимости для выполнения вычислений, связанных с диффурами
from scipy. integrate import odeint

# Инициализация приложения FastAPI
app = FastAPI()

# Вспомогательный класс F, выражающий полином третьей степени.
# Поля:
# a - коэффициент при x^3
# b - коэффициент при x^2
# c - коэффициент при x
# d - свободный коэффициент
# Методы:
# calc(x) - вычисляет значение полинома при заданном x
class F:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def calc(self, x):
        return self.a* (x**3) + self.b * (x**2) + self.c * x + self.d

# Функция, выражающая систему уравнений Практической Работы 2
# Принимает:
# u - массив исследуемых параметров, представляющих собой функции xi(t), i=[1..14]
# t - массив временных точек от 0 до 1
# с - словарь нормировочных множителей вида {"Название": <Значение множителя>, ...}
# f - словарь полиномов и возмущений вида {"Название": <Экземпляр полинома>, ...}
# Возвращает:
# массив вычисленных выражений dxi/dt, i=[1..14]
def du2_dt(u,t,c,f):

        # Извлекаем элементы, которые являются исследуемыми xi(t), из массива u
        [x1_t, x2_t, x3_t, x4_t, x5_t, x6_t, x7_t, x8_t, x9_t, x10_t, x11_t, x12_t, x13_t, x14_t] = u

        # Для каждого уравнения приводим запись формулы из задания. 
        # Для формулы:
        # dx1/dt = 1/X1norm * 
        # ( 
        #   f1(x1(t)) * f2(x2(t)) * f3(x3(t)) * f4(x4(t)) * f5(x5(t)) * f6(x6(t)) * f7(x7(t)) *
        #   f8(x8(t)) * f9(x9(t)) * f10(x10(t)) * f11(x11(t)) * f12(x12(t)) * f13(x13(t)) * f14(x14(t)) * 
        #   (z1(t) + z2(t) + z3(t))
        #   - z4(t) - z5(t)
        # )
        # программная запись примет вид:
        dx1_dt = (
            (1/c['normX1']) *
            (
                f['F1'].calc(x1_t) *
                f['F2'].calc(x2_t) *
                f['F3'].calc(x3_t) *
                f['F4'].calc(x4_t) *
                f['F5'].calc(x5_t) *
                f['F6'].calc(x6_t) *
                f['F7'].calc(x7_t) *
                f['F8'].calc(x8_t) *
                f['F9'].calc(x9_t) *
                f['F10'].calc(x10_t) *
                f['F11'].calc(x11_t) *
                f['F12'].calc(x12_t) *
                f['F13'].calc(x13_t) *
                f['F14'].calc(x14_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t)
                ) -
                f['Z4'].calc(t) -
                f['Z5'].calc(t)
            )
        )

        #Далее по аналогии со всеми уравнениями системы...
        dx2_dt = (
            (1/c['normX2']) *
            (
                f['F15'].calc(x1_t) *
                f['F16'].calc(x2_t) *
                f['F17'].calc(x3_t) *
                f['F18'].calc(x4_t) *
                f['F19'].calc(x5_t) *
                f['F20'].calc(x6_t) *
                f['F21'].calc(x7_t) *
                f['F22'].calc(x8_t) *
                f['F23'].calc(x9_t) *
                f['F24'].calc(x10_t) *
                f['F25'].calc(x11_t) *
                f['F26'].calc(x12_t) *
                f['F27'].calc(x13_t) *
                f['F28'].calc(x14_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t)
                ) -
                f['Z5'].calc(t)
            )
        )

        dx3_dt = (
            (1/c['normX3']) *
            (
                f['F29'].calc(x1_t) *
                f['F30'].calc(x2_t) *
                f['F31'].calc(x3_t) *
                f['F32'].calc(x4_t) *
                f['F33'].calc(x5_t) *
                f['F34'].calc(x6_t) *
                f['F35'].calc(x7_t) *
                f['F36'].calc(x8_t) *
                f['F37'].calc(x9_t) *
                f['F38'].calc(x10_t) *
                f['F39'].calc(x11_t) *
                f['F40'].calc(x12_t) *
                f['F41'].calc(x13_t) *
                f['F42'].calc(x14_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t)
                ) -
                f['Z5'].calc(t)
            )
        )

        dx4_dt = (
            (1/c['normX4']) *
            (
                f['F43'].calc(x1_t) *
                f['F44'].calc(x2_t) *
                f['F45'].calc(x3_t) *
                f['F46'].calc(x4_t) *
                f['F49'].calc(x7_t) *
                f['F50'].calc(x8_t) *
                f['F51'].calc(x9_t) *
                f['F52'].calc(x10_t) *
                f['F53'].calc(x11_t) *
                f['F54'].calc(x12_t) *
                f['F55'].calc(x13_t) *
                f['F56'].calc(x14_t) *
                f['Z5'].calc(t) -
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t) +
                    f['F47'].calc(x5_t) *
                    f['F48'].calc(x6_t)
                )
            )
        )

        dx5_dt = (
            (1/c['normX5']) *
            (
                f['F57'].calc(x4_t) *
                f['F58'].calc(x6_t) *
                f['F59'].calc(x9_t) *
                f['F60'].calc(x10_t) *
                f['F61'].calc(x13_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                ) -
                f['Z5'].calc(t)
            )
        )

        dx6_dt = (
            (1/c['normX6']) *
            (
                f['F62'].calc(x1_t) *
                f['F63'].calc(x2_t) *
                f['F64'].calc(x3_t) *
                f['F65'].calc(x4_t) *
                f['F66'].calc(x5_t) *
                f['F67'].calc(x6_t) *
                f['F68'].calc(x7_t) *
                f['F69'].calc(x8_t) *
                f['F70'].calc(x9_t) *
                f['F71'].calc(x10_t) *
                f['F72'].calc(x11_t) *
                f['F73'].calc(x12_t) *
                f['F74'].calc(x13_t) *
                f['F75'].calc(x14_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t)
                ) -
                f['Z5'].calc(t)
            )
        )

        dx7_dt = (
            (1/c['normX7']) *
            (
                f['F151'].calc(x2_t) *
                f['F152'].calc(x4_t) *
                f['F153'].calc(x14_t) -
                f['Z5'].calc(t)
            )
        )

        dx8_dt = (
            (1/c['normX8']) *
            (
                f['F76'].calc(x1_t) *
                f['F77'].calc(x2_t) *
                f['F78'].calc(x3_t) *
                f['F79'].calc(x4_t) *
                f['F80'].calc(x6_t) *
                f['F81'].calc(x9_t) *
                f['F82'].calc(x10_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t)
                ) -
                f['Z1'].calc(t) -
                f['Z2'].calc(t)
            )
        )

        dx9_dt = (
            (1/c['normX9']) *
            (
                f['F83'].calc(x1_t) *
                f['F84'].calc(x2_t) *
                f['F85'].calc(x3_t) *
                f['F86'].calc(x4_t) *
                f['F87'].calc(x5_t) *
                f['F88'].calc(x6_t) *
                f['F89'].calc(x7_t) *
                f['F90'].calc(x10_t) *
                f['F91'].calc(x11_t) *
                f['F92'].calc(x12_t) *
                f['F93'].calc(x13_t) *
                f['F94'].calc(x14_t) *
                (
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                ) -
                f['Z1'].calc(t) -
                f['Z2'].calc(t) -
                f['Z3'].calc(t)
            )
        )

        dx10_dt = (
            (1/c['normX10']) *
            (
                f['F95'].calc(x1_t) *
                f['F96'].calc(x2_t) *
                f['F97'].calc(x3_t) *
                f['F98'].calc(x4_t) *
                f['F99'].calc(x5_t) *
                f['F100'].calc(x6_t) *
                f['F101'].calc(x7_t) *
                f['F102'].calc(x8_t) *
                f['F103'].calc(x9_t) *
                f['F104'].calc(x10_t) *
                f['F105'].calc(x11_t) *
                f['F106'].calc(x12_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t)
                ) -
                f['F107'].calc(x13_t) *
                f['F108'].calc(x14_t) *
                f['Z3'].calc(t)
            )
        )

        dx11_dt = (
            (1/c['normX11']) *
            (
                f['F109'].calc(x1_t) *
                f['F110'].calc(x2_t) *
                f['F111'].calc(x3_t) *
                f['F112'].calc(x4_t) *
                f['F113'].calc(x8_t) *
                f['F114'].calc(x10_t) *
                f['F115'].calc(x12_t) *
                f['F116'].calc(x13_t) *
                f['F117'].calc(x14_t) -
                (
                    f['F118'].calc(x5_t) *
                    f['F119'].calc(x6_t)
                ) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                )
            )
        )

        dx12_dt = (
            (1/c['normX12']) *
            (
                f['F120'].calc(x1_t) *
                f['F121'].calc(x2_t) *
                f['F122'].calc(x3_t) *
                f['F123'].calc(x4_t) *
                f['F124'].calc(x5_t) *
                f['F125'].calc(x6_t) *
                f['F126'].calc(x7_t) *
                f['F127'].calc(x8_t) *
                f['F128'].calc(x9_t) *
                f['F129'].calc(x10_t) *
                f['F130'].calc(x11_t) *
                f['F131'].calc(x13_t) *
                f['F132'].calc(x14_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                ) -
                f['Z3'].calc(t)
            )
        )

        dx13_dt = (
            (1/c['normX13']) *
            (
                f['F133'].calc(x1_t) *
                f['F134'].calc(x2_t) *
                f['F135'].calc(x3_t) *
                f['F136'].calc(x4_t) *
                f['F137'].calc(x5_t) *
                f['F138'].calc(x6_t) *
                f['F139'].calc(x7_t) *
                f['F140'].calc(x8_t) *
                f['F141'].calc(x10_t) *
                f['F142'].calc(x11_t) *
                f['F143'].calc(x12_t) *
                f['F144'].calc(x13_t) *
                f['F145'].calc(x14_t) -
                f['F146'].calc(x9_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                )
            )
        )

        dx14_dt = (
            (1/c['normX14']) *
            (
                f['F147'].calc(x5_t) *
                f['F148'].calc(x7_t) *
                f['F149'].calc(x11_t) *
                f['F150'].calc(x13_t) -
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t)
                )
            )
        )

        # Возвращаем массив вычисленных dxi/dt, i=[1..14]
        return [dx1_dt, dx2_dt, dx3_dt, dx4_dt, dx5_dt, dx6_dt, dx7_dt, dx8_dt, dx9_dt, dx10_dt, dx11_dt, dx12_dt, dx13_dt,
        dx14_dt]

# Указывает маршрут /lab2 до web-интерфейса, созданного для ввода значений из Практической Работы 2
@app.get("/lab2")
async def main2():
    # Возвращаем в качестве ответа HTML-страницу, расположенную по адресу /home/AlexAranara/my_fastapi/web/lab2.html
    return FileResponse("venv/lab2.html")

# Указываем маршрут /api/lab2/count, по которому можно получить результаты вычислений для Практической работы 2
# Тело данного POST-запрос имеет вид {"x0": <Словарь начальный значений>, "c": <Словарь нормировочных множителей>, "f": <Словарь полиномов и возмущений>}
# Записи <Словарь начальных значений> имеют вид "Имя значения": <значение>
# Записи <Словарь нормировочных множителей> имеют вид "Имя нормировочного множителя": <значение>
# Записи <Словарь полиномов и возмущений> имеют вид "Имя полинома (возмущения)": {"a": <значение>, "b": <значение>, "c": <значение>, "d": <значение> }
@app.post("/api/lab2/count")
async def count2(data  = Body()):
    try:
        # Извлекаем пришедшие начальные значения в t0
        t0 = data['x0']

        # Устанавливаем значения t в диапазоне от 0 до 1 с шагом 0.05, при которых будут проводиться вычисления
        t_span = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]

        # Извлекаем пришедшие нормировочные множители в с
        c = data['c']

        # Извлекаем пришедшие коэффициенты полиномов и возмущений в uf
        uf = data['f']

        # На основе коэффициентов uf создадим удобные для работы экземпляры класса F, каждый из которых сопоставив с названием
        # Соберем пары "Имя полинома (возмущения)" : <Экземпляр класса F> в словарь f
        f = {}
        for key in uf:
            f[key] = F(uf[key]['a'],uf[key]['b'],uf[key]['c'],uf[key]['d'])

        # Передаем в качестве аргументов в функцию для решения системы диффуров odeint 
        # Параметры:
        # du_dt - ранее описанная функция, выражающая систему уравнений Практической Работы 1
        # list(t0.values()) - список начальных значений
        # t_span - массив временных точек от 0 до 1
        # args=(c,f,) - дополнительный аргументы в виде констант и полиномов
        # full_output=True - для определения, решена ли система или нет
        # Получаем на выходе:
        # usolution - Решение системы в виде массива, элементы которого являются массивами решений для указанной точке на временной прямой
        # d - лог вычислений
        usolution,d = odeint(du2_dt, list(t0.values()), t_span, args=(c,f,), full_output=True)

        # Если решение не определено для заданных параметров - вернуть ошибку
        if (d["message"] != "Integration successful."):
            raise HTTPException(status_code=500)

        # Преобразуем значения для отправки их в качетсве JSON-ответа
        solution = [None] * len(usolution)
        idx = 0
        for elem in usolution:
            solution[idx] = list(elem)
            idx = idx + 1

        return JSONResponse(content={"message": solution})
    except:
        # При некорректно предоставленных в теле запроса данных
        raise HTTPException(status_code=500)

