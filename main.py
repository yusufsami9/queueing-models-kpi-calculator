import tkinter as tk
from tkinter import ttk
from tkinter import StringVar
from tkinter import messagebox
import math

pn_dict = {}


def queuing_model_1(lam, mu):  # M/M/1:GD/∞/∞
    if lam >= mu:
        messagebox.showerror(
            "Input Error",
            "λ (arrival rate) should be less than µ (service rate). Please enter valid values."
        )
        return None

    p0 = (mu - lam) / mu
    Ls = lam / (mu - lam)
    Ws = Ls / lam
    Wq = Ws - (1 / mu)
    Lq = lam * Wq
    c_bar = lam / mu
    lam_eff = lam
    lam_lost = 0.0

    global pn_dict
    pn_dict = {}
    for n in range(1, 21):
        pn_dict[f"p_{n}"] = (lam / mu) ** n * p0

    return {
        "p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq,
        "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost
    }


def queuing_model_2(lam, mu, N):  # M/M/1:GD/N/∞
    if lam >= mu:
        messagebox.showerror(
            "Input Error",
            "λ (arrival rate) should be less than µ (service rate). Please enter valid values."
        )
        return None
    if N is None or N < 0:
        messagebox.showerror("Input Error", "N must be a non-negative integer.")
        return None

    rho = lam / mu
    if rho != 1:
        p0 = (1 - rho) / (1 - rho ** (N + 1))
        pN = (rho ** N) * p0
    else:
        p0 = 1 / (N + 1)
        pN = p0

    global pn_dict
    pn_dict = {}
    for n in range(1, 21):
        pn_dict[f"p_{n}"] = (rho ** n) * p0 if n <= N else 0.0

    Ls = sum(n * pn_dict[f"p_{n}"] for n in range(1, N + 1))
    Lq = Ls - (1 - p0)

    lam_lost = lam * pN
    lam_eff = lam - lam_lost

    if lam_eff <= 0:
        Ws = 0.0
        Wq = 0.0
    else:
        Ws = Ls / lam_eff
        Wq = Lq / lam_eff

    c_bar = Ls - Lq

    return {
        "p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq,
        "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost
    }


def queuing_model_3(lam, mu, c):  # M/M/c:GD/∞/∞
    if c is None or c <= 0:
        messagebox.showerror("Input Error", "c must be a positive integer.")
        return None
    if lam >= c * mu:
        messagebox.showerror(
            "Input Error",
            "λ (arrival rate) should be less than µ * c. Please enter valid values."
        )
        return None

    rho = lam / (c * mu)

    summation = sum((lam ** n) / (math.factorial(n) * (mu ** n)) for n in range(c))
    second_term = ((lam ** c) / (math.factorial(c) * (mu ** c))) * (1 / (1 - rho))
    p0 = 1 / (summation + second_term)

    # Standard Erlang-C based Lq for M/M/c
    Lq = (
        (p0 * (lam ** c) * rho) /
        (math.factorial(c) * (mu ** c) * ((1 - rho) ** 2))
    )

    Ls = Lq + (lam / mu)
    Wq = Lq / lam
    Ws = Wq + (1 / mu)
    c_bar = Ls - Lq

    global pn_dict
    pn_dict = {}
    for n in range(1, 21):
        if n <= c:
            pn_dict[f"p_{n}"] = ((lam ** n) / (math.factorial(n) * (mu ** n))) * p0
        else:
            pn_dict[f"p_{n}"] = (lam ** n / ((c ** (n - c)) * math.factorial(c) * (mu ** n))) * p0

    lam_eff = lam
    lam_lost = 0.0

    return {
        "p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq,
        "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost
    }


def queuing_model_4(lam, mu, c, N):  # M/M/c:GD/N/∞
    if c is None or c <= 0:
        messagebox.showerror("Input Error", "c must be a positive integer.")
        return None
    if N is None or N < c:
        messagebox.showerror("Input Error", "N must be an integer and N >= c.")
        return None
    if mu == 0:
        messagebox.showerror("Input Error", "µ must be > 0.")
        return None

    rho = lam / (c * mu)
    if rho == 1:
        messagebox.showerror("Input Error", "This implementation assumes ρ ≠ 1 for finite-capacity M/M/c/N.")
        return None

    summation_term = sum((lam ** n) / (math.factorial(n) * (mu ** n)) for n in range(c))

    second_term_numerator = (lam ** c) * (1 - (lam / (c * mu)) ** (N - c + 1))
    second_term_denominator = (math.factorial(c) * (mu ** c) * (1 - lam / (c * mu)))
    second_term = second_term_numerator / second_term_denominator

    p0 = 1 / (summation_term + second_term)

    term1 = (lam ** (c + 1))
    term2 = ((c - lam / mu) ** 2) * math.factorial(c - 1) * (mu ** (c + 1))
    term3 = 1 - (lam / (c * mu)) ** (N - c + 1)
    term4 = (N - c + 1) * (1 - (lam / (c * mu))) * (lam / (c * mu)) ** (N - c)

    Lq = (term1 / term2) * (term3 - term4) * p0

    global pn_dict
    pn_dict = {}
    for n in range(0, N + 1):
        if n <= c:
            pn = ((lam ** n) / (math.factorial(n) * (mu ** n))) * p0
        else:
            pn = (lam ** n / ((c ** (n - c)) * math.factorial(c) * (mu ** n))) * p0
        pn_dict[f"p_{n}"] = pn

    # Lost arrivals = λ * P(N)
    pN = pn_dict.get(f"p_{N}", 0.0)
    lam_lost = lam * pN
    lam_eff = lam - lam_lost

    Ls = Lq + (lam_eff / mu)

    if lam_eff > 0:
        Ws = Ls / lam_eff
        Wq = Lq / lam_eff
    else:
        Ws = 0.0
        Wq = 0.0

    c_bar = Ls - Lq

    return {
        "p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq,
        "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost
    }


def queuing_model_5(lam, mu):  # M/M/∞:GD/∞/∞ (infinite servers)
    if mu == 0:
        messagebox.showerror("Input Error", "µ must be > 0.")
        return None

    Lq = 0.0
    Wq = 0.0
    Ls = lam / mu
    Ws = 1 / mu
    c_bar = Ls
    lam_eff = lam
    lam_lost = 0.0
    p0 = math.exp(-lam / mu)

    global pn_dict
    pn_dict = {}
    for n in range(1, 21):
        pn_dict[f"p_{n}"] = ((lam ** n) / (math.factorial(n) * (mu ** n))) * p0

    return {
        "p0": p0, "Ls": Ls, "Ws": Ws, "Wq": Wq, "Lq": Lq,
        "c_bar": c_bar, "lam_eff": lam_eff, "lam_lost": lam_lost
    }


def hide_or_show(selected, c_label, c_entry, N_label, N_entry, models):
    if selected == models[0]:  # M/M/1 infinite
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid_forget()
        N_entry.grid_forget()
    elif selected == models[1]:  # M/M/1/N
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid(row=0, column=2, padx=20, pady=5, sticky="e")
        N_entry.grid(row=0, column=3, padx=10, pady=5, sticky="w")
    elif selected == models[2]:  # M/M/c infinite
        c_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        c_entry.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        N_label.grid_forget()
        N_entry.grid_forget()
    elif selected == models[3]:  # M/M/c/N
        c_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        c_entry.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        N_label.grid(row=0, column=2, padx=20, pady=5, sticky="e")
        N_entry.grid(row=0, column=3, padx=10, pady=5, sticky="w")
    elif selected == models[4]:  # M/M/infty
        c_label.grid_forget()
        c_entry.grid_forget()
        N_label.grid_forget()
        N_entry.grid_forget()


def show_warning():
    messagebox.showwarning(
        "Warning",
        "Error: Some variables are not filled or a non-numerical character has been entered."
    )


def update_values_table(values_table_frame, metrics):
    for widget in values_table_frame.winfo_children():
        widget.destroy()

    if metrics is None:
        return

    tree = ttk.Treeview(values_table_frame, columns=list(metrics.keys()), show="headings")

    for key in metrics.keys():
        tree.heading(key, text=key)
        tree.column(key, width=90, anchor=tk.CENTER)

    tree.insert("", tk.END, values=list(metrics.values()))
    tree.pack(fill=tk.BOTH, expand=True)


def update_p_table(p_table_frame, pn):
    for widget in p_table_frame.winfo_children():
        widget.destroy()

    tree = ttk.Treeview(p_table_frame, columns=("N", "P(N)"), show="headings")
    tree.heading("N", text="N")
    tree.heading("P(N)", text="P(N)")
    tree.column("N", width=120, anchor=tk.CENTER)
    tree.column("P(N)", width=220, anchor=tk.CENTER)

    for key, value in pn.items():
        try:
            tree.insert("", tk.END, values=(key, f"{float(value):.5f}"))
        except Exception:
            tree.insert("", tk.END, values=(key, str(value)))

    tree.pack(fill=tk.BOTH, expand=True)


def main():
    global pn_dict

    root = tk.Tk()
    root.title("Queuing Models")
    root.geometry("650x780")
    root.configure(bg="lightblue")

    models = [
        "M/M/1:GD/∞/∞",
        "M/M/1:GD/N/∞",
        "M/M/c:GD/∞/∞",
        "M/M/c:GD/N/∞",
        "M/M/∞:GD/∞/∞"
    ]

    selected_option = StringVar(value=models[0])

    tk.Label(
        root,
        text="Operation Research II – Queuing Models Project",
        bg="lightblue",
        wraplength=600,
        font=("Arial", 16, "bold")
    ).pack(pady=10)

    dropdown = ttk.Combobox(root, textvariable=selected_option, values=models, state="readonly")
    dropdown.pack(pady=15)

    label_font = ("Arial", 13, "bold")

    frame_lambda_mu = tk.Frame(root, bg="lightblue")
    frame_lambda_mu.pack(padx=20, pady=15)

    tk.Label(frame_lambda_mu, text="λ (arrival rate):", bg="lightblue", font=label_font).grid(row=0, column=0, padx=10, pady=5)
    lam_var = tk.StringVar()
    lam_entry = tk.Entry(frame_lambda_mu, textvariable=lam_var, width=20)
    lam_entry.grid(row=0, column=1, padx=10, pady=5)

    tk.Label(frame_lambda_mu, text="µ (service rate):", bg="lightblue", font=label_font).grid(row=0, column=2, padx=20, pady=5)
    mu_var = tk.StringVar()
    mu_entry = tk.Entry(frame_lambda_mu, textvariable=mu_var, width=20)
    mu_entry.grid(row=0, column=3, padx=10, pady=5)

    frame_c_N = tk.Frame(root, bg="lightblue")
    frame_c_N.pack(pady=5)

    c_label = tk.Label(frame_c_N, text="c (servers):", bg="lightblue", font=label_font)
    c_var = tk.StringVar()
    c_entry = tk.Entry(frame_c_N, textvariable=c_var, width=20)

    N_label = tk.Label(frame_c_N, text="N (capacity):", bg="lightblue", font=label_font)
    N_var = tk.StringVar()
    N_entry = tk.Entry(frame_c_N, textvariable=N_var, width=20)

    def on_selection_change(*args):
        hide_or_show(selected_option.get(), c_label, c_entry, N_label, N_entry, models)

    selected_option.trace_add("write", on_selection_change)
    hide_or_show(selected_option.get(), c_label, c_entry, N_label, N_entry, models)

    p_table_frame = ttk.Frame(root)
    p_table_frame.pack(fill=tk.BOTH, expand=False, padx=30, pady=20)

    values_table_frame = ttk.Frame(root)
    values_table_frame.pack(fill=tk.BOTH, expand=False, padx=30, pady=10)

    def calculate():
        try:
            lam = float(lam_var.get())
            mu = float(mu_var.get())
            c = int(c_var.get()) if c_var.get().strip() else None
            N = int(N_var.get()) if N_var.get().strip() else None

            selected = selected_option.get()

            metrics = None
            if selected == models[0]:
                metrics = queuing_model_1(lam, mu)
            elif selected == models[1]:
                metrics = queuing_model_2(lam, mu, N)
            elif selected == models[2]:
                metrics = queuing_model_3(lam, mu, c)
            elif selected == models[3]:
                metrics = queuing_model_4(lam, mu, c, N)
            elif selected == models[4]:
                metrics = queuing_model_5(lam, mu)

            update_p_table(p_table_frame, pn_dict)
            update_values_table(values_table_frame, metrics)

        except ValueError:
            show_warning()

    button = ttk.Button(root, text="Calculate", command=calculate)
    button.pack(pady=10)

    root.mainloop()


if __name__ == "__main__":
    main()
