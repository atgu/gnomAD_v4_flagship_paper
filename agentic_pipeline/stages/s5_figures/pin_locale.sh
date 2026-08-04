#!/usr/bin/env bash
# Pin the collation the committed figures were rendered with.
#
# Sourced, not executed, and it exits the caller on failure. Sets LC_ALL and then
# verifies that R honoured it, because a locale that is not generated on the
# machine is accepted silently by setlocale() and afterwards ignored — which
# would change a figure with no error reported anywhere.

PEPPER_LOCALE="${PEPPER_LOCALE:-en_US.UTF-8}"

# locale -a prints "en_US.utf8" where the environment variable is spelled
# "en_US.UTF-8", so compare with punctuation and case removed.
_norm() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | tr -d '._-'; }

_wanted="$(_norm "$PEPPER_LOCALE")"
_found=0
while IFS= read -r _l; do
  [ "$(_norm "$_l")" = "$_wanted" ] && { _found=1; break; }
done < <(locale -a 2>/dev/null)

if [ "$_found" -eq 0 ]; then
  cat >&2 <<EOF
ERROR: locale $PEPPER_LOCALE is not available on this machine.

The figures were rendered with its collation, and it decides where the "<2015"
category lands on Figure 6's panel A. Generate it, for instance with

    sudo locale-gen en_US.UTF-8 && sudo update-locale

or override the choice with PEPPER_LOCALE=<locale>, accepting that panel A may
be ordered differently from the committed reference.
EOF
  exit 1
fi

export LC_ALL="$PEPPER_LOCALE"

_actual="$(Rscript -e 'cat(Sys.getlocale("LC_COLLATE"))' 2>/dev/null)"
if [ "$(_norm "$_actual")" != "$_wanted" ]; then
  echo "ERROR: R reports LC_COLLATE=$_actual, expected $PEPPER_LOCALE." >&2
  echo "       Sorting would differ from the committed figures; refusing to run." >&2
  exit 1
fi

echo "Collation : $PEPPER_LOCALE (pinned)"

unset -f _norm
unset _wanted _found _l _actual
