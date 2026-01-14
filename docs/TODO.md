# add search PCRE like
- \E : stop interne (* avec au moins un non-gap à droite)
- \e : stop terminal (* au dernier résidu non-gap ; peut être suivi de -)
- \F : frameshift interne (! interne)
- \f : frameshift terminal (! terminal)
- \G : gap interne (- avant/au dernier non-gap)
- \g : gap terminal (- dans la traîne finale après le dernier non-gap)
Et en brut :
- \*, \!, \- : match caractère exact.
