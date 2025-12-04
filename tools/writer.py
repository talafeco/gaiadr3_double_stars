import os

def write_structured_output(filename_id, data_groups, folder="output"):
    """
    Kiírja a strukturált adatokat egy TXT fájlba, négy oszlopban: Név, Érték, ±, Hiba.
    Csoportosítva: SIMBAD, DR3, DR2, Caltech, Számított, Izokrón.
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    filename = os.path.join(folder, f"{filename_id}.txt")

    with open(filename, mode="w", encoding="utf-8") as f:
        f.write("Adatcsoport\tAdat neve\tÉrték\t±\tHiba\n")

        for group, items in data_groups.items():
            for name, (value, error) in items.items():
                pm = "±" if error not in ("", None, "--") else ""
                f.write(f"{group}\t{name}\t{value}\t{pm}\t{error}\n")

    print(f"💾 Adatok kiírva: {filename}")
