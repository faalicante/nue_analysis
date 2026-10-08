#!/usr/bin/env python3
"""Rename brick identifiers in file names.

Examples:
  ./rename_brick_files.py /path/to/b000222 22 222
  ./rename_brick_files.py /path/to/b000222 000022 000222 --apply
"""

from __future__ import annotations

import argparse
from pathlib import Path


def brick_token(value: str) -> str:
    """Return the canonical b-prefixed, six-digit brick token."""
    value = value.strip()
    if value.startswith("b"):
        value = value[1:]
    if not value.isdigit():
        raise argparse.ArgumentTypeError(f"brick id non numerico: {value!r}")
    return f"b{int(value):06d}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rinomina i file sostituendo un identificativo brick nel nome."
    )
    parser.add_argument("directory", type=Path, help="Cartella che contiene i file da rinominare")
    parser.add_argument(
        "from_brick",
        type=brick_token,
        help="Brick di partenza, es. 22, 022, 000022 o b000022",
    )
    parser.add_argument(
        "to_brick",
        type=brick_token,
        help="Brick di arrivo, es. 222, 000222 o b000222",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Esegue davvero la rinomina. Senza questa opzione mostra solo l'anteprima.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    directory = args.directory.expanduser()

    if not directory.is_dir():
        raise SystemExit(f"Cartella non trovata: {directory}")

    planned: list[tuple[Path, Path]] = []
    for source in sorted(directory.iterdir()):
        if not source.is_file() or args.from_brick not in source.name:
            continue
        target = source.with_name(source.name.replace(args.from_brick, args.to_brick))
        planned.append((source, target))

    if not planned:
        print(f"Nessun file contiene {args.from_brick} in {directory}")
        return 0

    collisions = [target for _, target in planned if target.exists()]
    if collisions:
        print("Rinomina bloccata: questi file di destinazione esistono gia':")
        for target in collisions:
            print(f"  {target}")
        return 1

    action = "Rinomino" if args.apply else "Anteprima"
    print(f"{action}: {args.from_brick} -> {args.to_brick}")
    for source, target in planned:
        print(f"  {source.name} -> {target.name}")

    if args.apply:
        for source, target in planned:
            source.rename(target)
        print(f"Fatto: rinominati {len(planned)} file.")
    else:
        print(f"Totale file da rinominare: {len(planned)}")
        print("Aggiungi --apply per eseguire davvero.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
