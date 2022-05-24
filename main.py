"""
def find_peak_bounds(data, point):
    delta = 10
    while True:
        if point + delta >= data.shape[0]:
            return data.shape[0] - point - 1
        a = [0] * 3
        a[0] = data[point - delta, 1]
        a[1] = (data[point, 1] - data[point - delta, 1]) / delta
        a[2] = (data[point + delta, 1] + data[point - delta, 1] - 2 * data[point, 1]) /\
               (2 * (delta ** 2))
        # f(x) = a[0] + a[1](x - point + delta) + a[2](x - point + delta)(x - point)
        for i in range(point - delta + 1, point + delta):
            if a[0] + a[1] * (i - point + delta) + a[2] * (i - point + delta) * (i - point) > 1.5 * data[i, 1]:
                return delta
        delta += 10
"""
from tkinter.filedialog import askopenfilename
import scipy.signal
import matplotlib.pyplot as plt
from tkinter import *
import numpy as np
from math import pi, exp, cos, sqrt, log
from unicodedata import lookup
from tkinter.messagebox import showinfo

root = Tk()

class PeakInfoDisplayer:
    curr = 0
    info = {}
    text_id = -1
    peaks = []
    properties = {}
    custom_base = 0
    base_id = -1
    edge_id = -1
    mid_id = -1
    new_edge = -1
    new_mid = -1
    curr_sym = IntVar()

    def __init__(self):
        self.info["area"] = []
        self.info["max_height"] = []
        self.info["width_on_half_height"] = []
        self.info["beta"] = []
        self.info["phi"] = []
        self.info["left"] = []
        self.info["mid"] = []
        self.info["right"] = []
        self.info["gamma"] = []
        self.info["sigma"] = []
        self.info["min_t"] = []
        self.info["lorenz_error"] = []
        self.info["gauss_error"] = []
        self.info["voigt_error"] = []
        self.info["lorenz_base"] = []
        self.info["gauss_base"] = []
        self.info["voigt_base"] = []
        self.info["lorenz_scale"] = []
        self.info["gauss_scale"] = []
        self.info["voigt_scale"] = []
        next_btn = Button(root, text='->', command=self.display_next)
        next_btn.place(x=1350, y=310)
        prev_btn = Button(root, text='<-', command=self.display_prev)
        prev_btn.place(x=1250, y=310)
        scan_new_peak_btn = Button(root, text="Выделить новый пик", command=self.scan_new_peak)
        scan_new_peak_btn.place(x=1057, y=80)
        add_new_peak_btn = Button(root, text="Добавить новый пик", command=self.add_new_peak)
        add_new_peak_btn.place(x=1057, y=110)
        delete_peak_btn = Button(root, text="Удалить этот пик", command=self.delete_peak)
        delete_peak_btn.place(x=1260, y=280)
        r1 = Radiobutton(root, text="Отразить слева направо", variable=self.curr_sym, value=0, command=self.change_symmetry)
        r2 = Radiobutton(root, text="Отразить справа налево", variable=self.curr_sym, value=1, command=self.change_symmetry)
        r1.place(x=1225, y=450)
        r2.place(x=1225, y=480)

    def delete_peak(self):
        for key in self.info.keys():
            self.info[key].pop(self.curr)
        if self.curr == len(self.info["mid"]):
            self.curr = 0
        if self.curr > 0:
            self.display_curr()

    def scan_new_peak(self):
        canvas.delete(self.base_id)
        canvas.bind("<Button-1>", self.set_edge)
        canvas.bind("<Button-3>", self.set_mid)

    def set_edge(self, event):
        if raw_data.shape[0] > 0:
            self.new_edge = int((event.x - left_margin) / (width - left_margin) * raw_data.shape[0])
            canvas.delete(self.edge_id)
            self.edge_id = canvas.create_line(event.x, 0, event.x, height - bot_margin, fill='blue', width=2)

    def set_mid(self, event):
        if raw_data.shape[0] > 0:
            self.new_mid = int((event.x - left_margin) / (width - left_margin) * raw_data.shape[0])
            canvas.delete(self.mid_id)
            self.mid_id = canvas.create_line(event.x, 0, event.x, height - bot_margin, fill='green', width=2)

    def add_new_peak(self):
        if self.new_edge != -1 and self.new_mid != -1:
            if self.new_edge > self.new_mid:
                right = self.new_edge
                left = max(0, 2 * self.new_mid - right)
            else:
                left = self.new_edge
                right = min(raw_data.shape[0], 2 * self.new_mid - left)
            self.process_peak(raw_data, left, self.new_mid, right)
            self.curr = len(self.info["mid"]) - 1
            self.display_curr()
        else:
            showinfo(title="Информация", message="Сначала выберите границу и центр нового пика")

    def display_curr(self):
        text = ""
        if len(self.info["mid"]) > 0:
            text += "Обрабатывается пик в точке " + str(round(raw_data[self.info["mid"][self.curr], 0], 3))
            text += "\nПлощадь = " + self.info["area"][self.curr]
            text += "\nШирина на половине высоты = " + self.info["width_on_half_height"][self.curr]
            text += "\nМаксимальная высота = " + self.info["max_height"][self.curr]
            text += "\n" + lookup("GREEK SMALL LETTER BETA") + " = " + self.info["beta"][self.curr]
            text += "\n" + lookup("GREEK SMALL LETTER PHI") + " = " + self.info["phi"][self.curr]
            text += "\nПриближение функцией Гаусса:\n" + str(round(self.info["gauss_base"][self.curr], 3))
            text += " + " + str(round(self.info["gauss_scale"][self.curr] / self.info["sigma"][self.curr], 3)) + "/"
            text += "sqrt(2" + lookup("GREEK SMALL LETTER PI") + ") *\n exp(-1/2 * ((x - "
            text += str(round(raw_data[self.info["mid"][self.curr], 0], 3))
            text += ") / " + str(round(self.info["sigma"][self.curr], 3)) + ")^2)\nСреднеквадратичная ошибка: "
            text += str(round(self.info["gauss_error"][self.curr], 3)) + "\n"
            text += "Приближение функцией Лоренца:\n" + str(round(self.info["lorenz_base"][self.curr], 3)) + " + "
            text += str(round(self.info["lorenz_scale"][self.curr] / self.info["gamma"][self.curr], 3))
            text += " / (" + lookup("GREEK SMALL LETTER PI") + " * \n(1 + ((x - "
            text += str(round(raw_data[self.info["mid"][self.curr], 0], 3))
            text += ") / " + str(round(self.info["gamma"][self.curr], 3)) + ")^2))\n"
            text += "Среднеквадратичная ошибка: " + str(round(self.info["lorenz_error"][self.curr], 3)) + "\n"
            text += "Приближение псевдо-Фойгтом:\n" + str(round(self.info["voigt_base"][self.curr], 3)) + " + "
            text += str(round(self.info["voigt_scale"][self.curr], 3)) + " * \n(" +\
                    str(round(self.info["min_t"][self.curr], 3))
            text += "G(" + str(round(raw_data[self.info["mid"][self.curr], 0], 3)) + ", "
            text += str(round(self.info["sigma"][self.curr], 3)) + ") + " + str(round(1 - self.info["min_t"][self.curr],
                                                                                      3))
            text += "L(" + str(round(raw_data[self.info["mid"][self.curr], 0], 3)) + ", "
            text += str(round(self.info["gamma"][self.curr], 3)) + "))\nСредневадратичная ошибка: "
            text += str(round(self.info["voigt_error"][self.curr], 3)) + "\n"
        text_canvas.delete(self.text_id)
        self.text_id = text_canvas.create_text(260, 150, text=text, font='TkMenuFont', fill='black', justify='right')

    def change_symmetry(self):
        left = self.info["left"][self.curr]
        mid = self.info["mid"][self.curr]
        right = self.info["right"][self.curr]
        base = min(raw_data[left:right, 1])
        self.info["gamma"][self.curr], self.info["lorenz_scale"][self.curr], \
        self.info["lorenz_base"][self.curr], self.info["lorenz_error"][self.curr] = \
            lorenz_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.info["sigma"][self.curr], self.info["gauss_scale"][self.curr], \
        self.info["gauss_base"][self.curr], self.info["gauss_error"][self.curr] = \
            gauss_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.info["min_t"][self.curr], self.info["voigt_scale"][self.curr], \
        self.info["voigt_base"][self.curr], self.info["voigt_error"][self.curr] = \
            voigt_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.display_curr()

    def display_next(self):
        self.curr = (self.curr + 1) % len(self.info["mid"])
        self.display_curr()

    def display_prev(self):
        self.curr = (self.curr - 1) % len(self.info["mid"])
        self.display_curr()

    def detect_peaks(self):
        self.curr = 0
        for key in self.info.keys():
            self.info[key].clear()
        self.peaks, self.properties = scipy.signal.find_peaks(raw_data[:, 1],
                                                              prominence=0.7 * min(
                                                                  raw_data[0:raw_data.shape[0] // 4, 1]))
        print(self.peaks, self.properties)
        for i in range(len(self.peaks)):
            if i == 0 or self.peaks[i] - self.peaks[i - 1] > 50:
                self.process_peak(raw_data, self.properties["left_bases"][i], self.peaks[i],
                                  min(raw_data.shape[0] - 1, 2 * self.peaks[i] - self.properties["left_bases"][i]))
        self.display_curr()

    def process_peak(self, data, left, mid, right):
        self.info["left"].append(left)
        self.info["mid"].append(mid)
        self.info["right"].append(right)
        base = data[left, 1]
        area = find_area(data, left, right, base)
        self.info["area"].append(str(round(area, 3)))
        zoom = data[left:right, :]
        peak = [mid - left]
        half_height_point = 0
        for i in range(left, mid):
            if data[i, 1] - base > (data[mid, 1] - base) / 2:
                half_height_point = i
                break
        width_on_half_height = 2 * (data[mid, 0] - data[half_height_point, 0])
        self.info["width_on_half_height"].append(str(round(width_on_half_height, 3)))
        max_height = max(zoom[:, 1])
        self.info["max_height"].append(str(round(max_height, 3)))
        self.info["beta"].append(str(round(area / max_height, 3)))
        self.info["phi"].append(str(round(width_on_half_height * max_height / area, 3)))
        g, l_s, l_b, l_e = lorenz_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        s, g_s, g_b, g_e = gauss_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        t, v_s, v_b, v_e = voigt_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.info["gamma"].append(g)
        self.info["lorenz_scale"].append(l_s)
        self.info["lorenz_base"].append(l_b)
        self.info["lorenz_error"].append(l_e)
        self.info["sigma"].append(s)
        self.info["gauss_scale"].append(g_s)
        self.info["gauss_base"].append(g_b)
        self.info["gauss_error"].append(g_e)
        self.info["min_t"].append(t)
        self.info["voigt_scale"].append(v_s)
        self.info["voigt_base"].append(v_b)
        self.info["voigt_error"].append(v_e)

    def scan_base(self):
        """
        base = float(entry.get())
        data = raw_data
        left = self.info["left"][self.curr]
        if data[left, 1] > base:
            while data[left, 1] > base:
                left -= 1
        else:
            while data[left, 1] < base:
                left += 1
        mid = self.info["mid"][self.curr]
        right = min(data.shape[0], 2 * mid - left)
        area = find_area(data, left, right, base)
        self.info["area"][self.curr] = (str(round(area, 3)))
        zoom = data[left:right, :]
        max_height = max(zoom[:, 1])
        x = 0
        while data[x, 1] < (max_height + base) / 2:
            x += 1
        width_on_half_height = 2 * (mid - x)
        self.info["width_on_half_height"][self.curr] = (str(round(width_on_half_height, 3)))
        self.info["beta"][self.curr] = (str(round(area / max_height, 3)))
        self.info["phi"][self.curr] = (str(round(width_on_half_height * max_height / area, 3)))
        """
        canvas.delete(self.edge_id)
        canvas.delete(self.mid_id)
        self.new_edge = -1
        self.new_mid = -1
        canvas.bind("<Button-1>", self.process_with_custom_base)
        text_canvas.create_text(1200, 190, text="Выберите фон мышкой")

    def process_with_custom_base(self, event):
        base = (height - bot_margin - event.y) * max(raw_data[:, 1]) / (height - bot_margin - top_margin)
        print(event.y, (height - bot_margin - event.y) * max(raw_data[:, 1]) / (height - bot_margin - top_margin))
        canvas.delete(self.base_id)
        self.base_id = canvas.create_line(0, event.y, width, event.y, fill='blue', width=2)
        data = raw_data
        left = self.info["left"][self.curr]
        if data[left, 1] > base:
            while data[left, 1] > base and left > 0:
                left -= 1
        else:
            while data[left, 1] < base and left < data.shape[0] - 1:
                left += 1
        mid = self.info["mid"][self.curr]
        right = min(data.shape[0], 2 * mid - left)
        print(data[right, 0], data[right, 1])
        if data[right, 1] > base:
            while data[right, 1] > base and right < data.shape[0] - 1:
                right += 1
        else:
            while data[right, 1] < base and right > 0:
                right -= 1
        print(left, data[left, 0], data[left, 1], data[right, 0], data[left, 0], data[mid, 0], data[mid, 1], mid,
              data.shape, base)
        area = find_area(data, left, right, base)
        self.info["area"][self.curr] = (str(round(area, 3)))
        zoom = data[left:right, :]
        max_height = max(zoom[:, 1])
        x = 0
        while data[x, 1] < (max_height + base) / 2:
            x += 1
        width_on_half_height = 2 * (mid - x) * (data[1, 0] - data[0, 0])
        self.info["width_on_half_height"][self.curr] = (str(round(width_on_half_height, 3)))
        self.info["beta"][self.curr] = (str(round(area / max_height, 3)))
        self.info["phi"][self.curr] = (str(round(width_on_half_height * max_height / area, 3)))
        self.info["gamma"][self.curr], self.info["lorenz_scale"][self.curr], \
        self.info["lorenz_base"][self.curr], self.info["lorenz_error"][self.curr] = \
            lorenz_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.info["sigma"][self.curr], self.info["gauss_scale"][self.curr], \
        self.info["gauss_base"][self.curr], self.info["gauss_error"][self.curr] = \
            gauss_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.info["min_t"][self.curr], self.info["voigt_scale"][self.curr], \
        self.info["voigt_base"][self.curr], self.info["voigt_error"][self.curr] = \
            voigt_approx(raw_data, left, mid, right, self.curr_sym.get(), base)
        self.display_curr()


# https://habr.com/ru/post/151157/
def gaussian_blur(data):
    a = [0] * 3
    sigma = 2
    sigma_inv_4 = sigma * sigma
    sigma_inv_4 = 1 / (sigma_inv_4 * sigma_inv_4)
    coef_a = sigma_inv_4 * (sigma * (sigma * (sigma * 1.1442707 + 0.0130625) - 0.7500910) + 0.2546730)
    coef_w = sigma_inv_4 * (sigma * (sigma * (sigma * 1.3642870 + 0.0088755) - 0.3255340) + 0.3016210)
    coef_b = sigma_inv_4 * (sigma * (sigma * (sigma * 1.2397166 - 0.0001644) - 0.6363580) - 0.0536068)
    z0_abs = exp(coef_a)
    z0_real = z0_abs * cos(coef_w)
    z2 = exp(coef_b)
    z0_abs_2 = z0_abs * z0_abs
    a[2] = 1 / (z2 * z0_abs_2)
    a[0] = (z0_abs_2 + 2 * z0_real * z2) * a[2]
    a[1] = -(2 * z0_real + z2) * a[2]
    b0 = 1 - a[0] - a[1] - a[2]
    res = data.copy()
    for i in range(3, data.shape[0]):
        res[i, 1] = b0 * data[i, 1] + a[0] * res[i - 1, 1] + a[1] * res[i - 2, 1] + a[2] * res[i - 3, 1]
    for i in range(data.shape[0] - 3, -1):
        res[i, 1] = b0 * data[i, 1] + a[0] * res[i - 1, 1] + a[1] * res[i - 2, 1] + a[2] * res[i - 3, 1]
    return res


def plot_graph_plt(data):
    x_data = data[:, 0]
    y_data = data[:, 1]
    plt.style.use('bmh')
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title('Зависимость интенсивности от угла')
    ax.set_ylabel('Интенсивность в у.е.')
    ax.set_xlabel('Угол(2 Тета) в градусах')
    ax.plot(x_data, y_data)
    ax.set_xlim(xmin=x_data[0], xmax=x_data[-1])
    plt.show()


def plot_graph_tkinter(data, clear=True):
    if clear:
        canvas.delete("all")
    max = data[0, 1]
    min = data[0, 1]
    for i in range(data.shape[0]):
        if data[i, 1] > max:
            max = data[i, 1]
        if data[i, 1] < min:
            min = data[i, 1]
    canvas.create_line(left_margin, 0, left_margin, height - bot_margin, width=2, fill='black')
    canvas.create_line(left_margin, height - bot_margin, width, height - bot_margin, width=2, fill='black')
    for i in range(0, round(max), round(max // 10)):
        y = top_margin + (max - i) * (height - bot_margin - top_margin) / max
        canvas.create_text(15, y, text=i, font='TkMenuFont', fill='black', justify='right')
        canvas.create_line(left_margin, y, width, y, width=1, fill='black')
    for i in range(11):
        canvas.create_text(left_margin + (width - left_margin) * i / 10, height - bot_margin // 2,
                           text=round(data[0, 0] + (data[-1, 0] - data[0, 0]) * i / 10),
                           font='TkMenuFont', fill='black', justify='right')
        canvas.create_line(left_margin + (width - left_margin) * i / 10, 0,
                           left_margin + (width - left_margin) * i / 10, height - bot_margin,
                           width=1, fill='black')
    for i in range(data.shape[0] - 1):
        canvas.create_line(left_margin + i * (width - left_margin) / (data.shape[0] - 1),
                           height - bot_margin - data[i, 1] * (height - top_margin - bot_margin) / max,
                           left_margin + (i + 1) * (width - left_margin) / (data.shape[0] - 1),
                           height - bot_margin - data[i + 1, 1] * (height - top_margin - bot_margin) / max,
                           width=2, fill='red')


def find_area(data, left, right, base):
    if left + 1 == data.shape[0]:
        return 0
    area = 0
    delta = data[left + 1, 0] - data[left, 0]
    for i in range(left, right):
        area += data[i, 1] - base
    return delta * area


# copy left to right if sym == 1, else copy right to left
def lorenz_approx(data, left, mid, right, sym, inp_base):
    min_error = 100000000
    min_gamma = 0
    min_scale = 0
    min_base = 0
    zoom = data[left:right, :]
    if sym == 1:
        for j in range(zoom.shape[0] // 2, zoom.shape[0]):
            zoom[j, 1] = zoom[zoom.shape[0] - j, 1]
    else:
        for j in range(0, zoom.shape[0] // 2):
            zoom[j, 1] = zoom[zoom.shape[0] - 1 - j, 1]
    for k in range(5, 10):
        base = k * inp_base / 10
        half_height_point = 0
        for i in range(zoom.shape[0]):
            if zoom[i, 1] - base > (data[mid, 1] - base) / 2:
                half_height_point = i
                break
        gamma = data[mid, 0] - data[left + half_height_point, 0]
        approx = np.ndarray(zoom.shape)
        for i in range(approx.shape[0]):
            approx[i, 0] = zoom[i, 0]
            approx[i, 1] = base + gamma / (pi * (gamma * gamma + (data[i + left, 0] - data[mid, 0])
                                                 * (data[i + left, 0] - data[mid, 0])))
        scale = (data[mid, 1] - base) / (approx[approx.shape[0] // 2, 1] - base)
        for i in range(approx.shape[0]):
            approx[i, 1] -= base
            approx[i, 1] *= scale
            approx[i, 1] += base
        '''
        for i in range(approx.shape[0] - 1):
            canvas.create_line(left_margin + (left + i) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - approx[i, 1] * (height - top_margin - bot_margin) / approx[approx.shape[0] // 2, 1],
                               left_margin + (left + (i + 1)) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - approx[i + 1, 1] * (height - top_margin - bot_margin) / approx[approx.shape[0] // 2, 1],
                               width=2, fill='red')
        for i in range(zoom.shape[0] - 1):
            canvas.create_line(left_margin + (left + i) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - zoom[i, 1] * (height - top_margin - bot_margin) / approx[approx.shape[0] // 2, 1],
                               left_margin + (left + (i + 1)) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - zoom[i + 1, 1] * (height - top_margin - bot_margin) / approx[approx.shape[0] // 2, 1],
                               width=2, fill='green')
        '''
        error = 0
        for i in range(approx.shape[0]):
            error += (approx[i, 1] - zoom[i, 1]) * (approx[i, 1] - zoom[i, 1])
        if error < min_error:
            min_error = error
            min_gamma = gamma
            min_scale = scale
            min_base = base
    min_error /= zoom.shape[0]
    return min_gamma, min_scale, min_base, min_error


def gauss_approx(data, left, mid, right, sym, inp_base):
    min_error = 100000000
    min_sigma = 0
    min_scale = 0
    min_base = 0
    zoom = data[left:right, :]
    if sym == 1:
        for j in range(zoom.shape[0] // 2, zoom.shape[0]):
            zoom[j, 1] = zoom[zoom.shape[0] - j, 1]
    else:
        for j in range(0, zoom.shape[0] // 2):
            zoom[j, 1] = zoom[zoom.shape[0] - 1 - j, 1]
    for k in range(5, 10):
        base = k * inp_base / 10
        half_height_point = 0
        for i in range(zoom.shape[0]):
            if zoom[i, 1] - base > (data[mid, 1] - base) / 2:
                half_height_point = i
                break
        sigma = (data[mid, 0] - data[left + half_height_point, 0]) / sqrt(2 * log(2))
        approx = np.ndarray(zoom.shape)
        for i in range(approx.shape[0]):
            approx[i, 0] = zoom[i, 0]
            approx[i, 1] = base + exp(-(zoom[i, 0] - data[mid, 0]) * (zoom[i, 0] - data[mid, 0]) / (2 * sigma * sigma))\
                           / (sigma * sqrt(2 * pi))
        scale = (data[mid, 1] - base) / (approx[approx.shape[0] // 2, 1] - base)
        for i in range(approx.shape[0]):
            approx[i, 1] -= base
            approx[i, 1] *= scale
            approx[i, 1] += base
        '''
        for i in range(approx.shape[0] - 1):
            canvas.create_line(left_margin + (left + i) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - approx[i, 1] * (height - top_margin - bot_margin) / approx[
                                   approx.shape[0] // 2, 1],
                               left_margin + (left + (i + 1)) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - approx[i + 1, 1] * (height - top_margin - bot_margin) / approx[
                                   approx.shape[0] // 2, 1],
                               width=2, fill='red')
        for i in range(zoom.shape[0] - 1):
            canvas.create_line(left_margin + (left + i) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - zoom[i, 1] * (height - top_margin - bot_margin) / approx[
                                   approx.shape[0] // 2, 1],
                               left_margin + (left + (i + 1)) * width / (raw_data.shape[0] - 1),
                               height - bot_margin - zoom[i + 1, 1] * (height - top_margin - bot_margin) / approx[
                                   approx.shape[0] // 2, 1],
                               width=2, fill='green')
        '''
        error = 0
        for i in range(approx.shape[0]):
            error += (approx[i, 1] - zoom[i, 1]) * (approx[i, 1] - zoom[i, 1])
        if error < min_error:
            min_error = error
            min_sigma = sigma
            min_scale = scale
            min_base = base
    min_error /= zoom.shape[0]
    return min_sigma, min_scale, min_base, min_error


def voigt_approx(data, left, mid, right, sym, inp_base):
    res_scale = 0
    min_base = 0
    res_error = 100000000
    res_t = 0
    zoom = data[left:right, :]
    if sym == 1:
        for j in range(zoom.shape[0] // 2, zoom.shape[0]):
            zoom[j, 1] = zoom[zoom.shape[0] - j, 1]
    else:
        for j in range(0, zoom.shape[0] // 2):
            zoom[j, 1] = zoom[zoom.shape[0] - 1 - j, 1]
    for k in range(5, 10):
        base = inp_base * k / 10
        half_height_point = 0
        for i in range(zoom.shape[0]):
            if zoom[i, 1] - base > (data[mid, 1] - base) / 2:
                half_height_point = i
                break
        sigma = (data[mid, 0] - data[left + half_height_point, 0]) / sqrt(2 * log(2))
        gamma = data[mid, 0] - data[left + half_height_point, 0]
        approx = np.ndarray(zoom.shape)
        min_t = 0
        min_error = 10000000
        min_scale = 0
        for j in range(11):
            t = j / 10
            for i in range(approx.shape[0]):
                approx[i, 0] = zoom[i, 0]
                approx[i, 1] = base + t * exp(-(zoom[i, 0] - data[mid, 0]) * (zoom[i, 0] - data[mid, 0]) /
                                              (2 * sigma * sigma)) / (sigma * sqrt(2 * pi)) + (1 - t) * gamma / (
                                       pi * (gamma * gamma +
                                             (data[i + left, 0] - data[mid, 0]) * (
                                                     data[i + left, 0] - data[mid, 0])))
            scale = (data[mid, 1] - base) / (approx[approx.shape[0] // 2, 1] - base)
            for i in range(approx.shape[0]):
                approx[i, 1] -= base
                approx[i, 1] *= scale
                approx[i, 1] += base
            error = 0
            for i in range(approx.shape[0]):
                error += (approx[i, 1] - zoom[i, 1]) * (approx[i, 1] - zoom[i, 1])
            if error < min_error:
                min_t = t
                min_error = error
                min_scale = scale
        if min_error < res_error:
            res_error = min_error
            res_scale = min_scale
            res_t = min_t
            min_base = base
    res_error /= zoom.shape[0]
    return res_t, res_scale, min_base, res_error


# https://python-scripts.com/tkinter
def get_data():
    global raw_data
    filepath = askopenfilename(filetypes=[("All Files", "*.*")])
    if not filepath:
        return
    f = open(filepath, 'r')
    '''
    line = f.readline()
    delimeter = ' ' * (line.count(' ') - 2)
    f.close()
    raw_data = np.loadtxt(filepath, delimiter=delimeter)
    plot_graph_tkinter(raw_data)
    '''
    lines = f.readlines()
    # raw_data.reshape((len(lines), 2))
    raw_data = np.ndarray((len(lines), 2))
    for i in range(len(lines)):
        if lines[i][0].isdigit():
            x, y = lines[i].split()
            raw_data[i, 0] = float(x)
            raw_data[i, 1] = float(y)
    # raw_data = gaussian_blur(raw_data)
    plot_graph_tkinter(raw_data)
    pid.detect_peaks()


def set_bounds(event):
    if raw_data.shape[0] > 0:
        global zoom_right, zoom_left, zoom_mid, left_id, right_id
        position = "(x={}, y={})".format(event.x, event.y)
        real_x = (event.x - left_margin) / (width - left_margin) * raw_data.shape[0]
        if zoom_mid == -1:
            zoom_left = real_x
            zoom_right = real_x
        elif real_x < zoom_mid:
            zoom_left = real_x
            zoom_right = min(raw_data.shape[0] - 1, 2 * zoom_mid - real_x)
        else:
            zoom_right = real_x
            zoom_left = max(0, 2 * zoom_mid - real_x)
        canvas.delete(left_id)
        left_id = canvas.create_line(event.x, 0, event.x, height - bot_margin, fill='blue', width=2)
        print(event.type, "event", position, real_x, raw_data[round(real_x), 0])


def set_mid(event):
    if raw_data.shape[0] > 0:
        global zoom_mid, zoom_right, zoom_left, mid_id
        position = "(x={}, y={})".format(event.x, event.y)
        real_x = (event.x - left_margin) / (width - left_margin) * raw_data.shape[0]
        if zoom_left == -1:
            zoom_mid = real_x

        canvas.delete(mid_id)
        mid_id = canvas.create_line(event.x, 0, event.x, height - bot_margin, fill='green', width=2)
        print(event.type, "event2", position)


def blur_data():
    global raw_data
    raw_data = gaussian_blur(raw_data)
    plot_graph_tkinter(raw_data)
    pid.detect_peaks()


if __name__ == "__main__":
    root.title('Анализ дифрактограмм')
    root.geometry('1420x620')

    width = 1020
    height = 620
    text_canvas = Canvas(root, width=1420 - width - 20, height=height)
    canvas = Canvas(root, width=width + 20, height=height)
    # canvas.bind("<Button-1>", set_bounds)
    # canvas.bind("<Button-3>", set_mid)

    left_margin = 30
    top_margin = 50
    bot_margin = 20

    #zoom_left = -1
    #zoom_mid = -1
   # zoom_right = -1
    #left_id = -1
   # mid_id = -1
   # right_id = -1
    # fr_buttons = Frame(root)
    # btn_open = Button(fr_buttons, text="Открыть", command=get_data())
    # btn_open.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
    # fr_buttons.grid(row=0, column=0, sticky="ns")
    pid = PeakInfoDisplayer()
    open_file_btn = Button(root, text='Открыть файл', command=get_data)
    open_file_btn.place(x=1075, y=20)
    blur_btn = Button(root, text='Применить Гауссово размытие', command=blur_data)
    blur_btn.place(x=1025, y=50)

    # btn.pack(side=RIGHT, padx=100, pady=200)
    raw_data = np.ndarray((0, 0))
    get_data()
    # blurred_data = gaussian_blur(raw_data)
    # X_data = blurred_data[:, 0]
    # Y_data = blurred_data[:, 1]
    # lorenz_approx(blurred_data, 2450, 2575, 1, 0.9 * blurred_data[2450, 1])
    canvas.place(x=0, y=0)
    # entry = Entry(root)
    # entry.place(x=1250, y=150)
    text_canvas.place(x=width + 20, y=0)
    change_base_btn = Button(root, text="Пересчитать со своим фоном", command=pid.scan_base)
    change_base_btn.place(x=1223, y=400)
    text_canvas.create_text(270, 370,
                            text="Чтобы пересчитать со своим фоном, нажмите на эту кнопку,"
                                 " а затем выберите фон левой кнопкой мыши",
                            font='TkMenuFont', fill='black', justify='center', width=200)
    text_canvas.create_text(77, 205,
                            text="Чтобы добавить новый пик, нажмите на кнопку \"Выделить новый пик\","
                                 "затем выберите край пика левой кнопкой мыши и центр правой кнопкой мыши, "
                                 "затем нажмите кнопку \"Добавить новый пик\"",
                            font='TkMenuFont', fill='black', justify='center', width=150)
    text_canvas.create_text(270, 440, text="Симметрия:",
                            font='TkMenuFont', fill='black', justify='center', width=150)
    # plot_graph_plt(blurred_data)
    root.mainloop()
