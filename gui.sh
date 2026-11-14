#!/usr/bin/env bash
# Atalho para abrir a interface gráfica do MiniSPICE.
# Uso: ./gui.sh            (do diretório do projeto)
#      ./gui.sh --port N   (porta alternativa)
cd "$(dirname "$0")"

if [ ! -x .venv/bin/minispice ] || ! .venv/bin/python -c 'import minispice' >/dev/null 2>&1; then
    echo "configurando o ambiente..."
    python3 -m venv .venv || exit 1
    .venv/bin/python -m pip install -e ".[dev]" || exit 1
fi

exec .venv/bin/minispice --gui "$@"
