# MRI_IR_Calculator_MatLab

## Documentation README

### Appendix: Equation syntax

For the presentation, converting equations back and forth between LibreOffice changed their syntax.
Example: Original simpler syntax like `M_0` would become `M rSub{ size 8{0}}`, and `newline` -> `{} #`.
Also, words and names weren't always automatically italicized so it had to be specified in some places.
This equation syntax that worked in the end, across platforms:
```
........................................................................................
Slide 03: FLAIR – Fluid Attenuated Inv. Recovery
........................................................................................
 size 12{ alignc{ stack{
E rSub{size8{t}} `≝` e rSup{size8{-`{ t over ital"T1" }}} 
~;~ M rSub{size8{t}} `=` M rSub{size8{0}} 
`-` left( M rSub{size8{0}} `-` M rSub{size8{1}} right) E rSub{size8{t}}
 {} # {} # 
bold{S} `≝` { M rSub{size8{ital"TI"}} } over { M rSub{size8{0}} } 
`=` 1 `-` left("1"`-`"-1"right) E rSub{size8{ital"TI"}} 
bold{ `=` underline{ 1 `-` 2E rSub{size8{ital"TI"}} } }
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}}`=`0 ~drarrow~ E rSub{size8{ital"TI"}} 
`=` e rSup{size8{ -`{ ital"TI" over ital"T1" } }}`=`{ 1 over 2 }
 {} # {} # drarrow~
bold{ ital"TI" `=` underline{ ital"T1" `·` "ln"left( 2 right) } }
 {} # {} # {} # alignr size8{ital"Q.E.D."}
} } }

........................................................................................
Slide 11: IR-Calc (0th: In Real Time)
........................................................................................
 size 12{ alignc{ stack{
E rSub{size8{t}} `≝` e rSup{size8{-`{ t over ital"T1" }}} 
~;~ M rSub{size8{t}} `=` M rSub{size8{0}} 
`-` left( M rSub{size8{0}} `-` M rSub{size8{1}} right) E rSub{size8{t}}
 {} # {} # 
{ M rSub{size8{1}} } over { M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left(ital"TR" - ital"TI"right) }} 
`=` 1`-`{ {E rSub{size8{ital"TR"}}} over {E rSub{size8{ital"TI"}}} }
 {} # {} # 
bold{S} `≝` { M rSub{size8{ital"TI"}} } over { M rSub{size8{0}} } `=` 1`-`
left( 1+{ { M rSub{size8{1}} } over { M rSub{size8{0}} } } right)E rSub{size8{ital"TI"}}
bold{ `=` {underline { left( 1+ital"E" rSub{size8{ital"TR"}} right) 
	`-`left("1+1"right)E rSub{size8{ital"TI"}}} } }
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}}`=`0 ~drarrow~ E rSub{size8{ital"TIn"}} 
`=` e rSup{size8{ -`{ ital"TIn" over ital "T1" } }} 
`=` { { 1+ital"E" rSub{size8{ital"TR"}} } over 2 }
 {} # {} # drarrow~ 
bold{ ital"TI" rSub{size8{n}} `=` underline{ 
ital"T1" left[`"ln"left( 2 right)
	`-`"ln"left( 1+ital"E" rSub{size8{ital"TR"}} right)` right ]`;~
E rSub{size8{ital"TR"}} `=` e rSup{size8{ -`{ ital"TR" over ital"T1" } }} } }
 {} # {} # {} # {} # alignr size8{"Q.E.D."}
} } }

........................................................................................
Slide 12: IR-Calc (1st: Effective TR)
........................................................................................
 size 12{ alignc{ stack{
E rSub{size8{t}} `≝` e rSup{size8{-`{ t over ital"T1" }}}
 {} # {} # 
{ M rSub{size8{1}} }over{ M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left(ital"TR"- ital"TI" right) }} 
`=` 1`-`{ E rSub{size8{ital"TR"}} }over{ E rSub{size8{ital"TI"}} } 
~;~ left( ital"TR" rightarrow ital"TR" rSub{size8{ital"eff"}} right)
 {} # {} # 
bold{S rSub{size8{0}}} `≝` { M rSub{size8{ital"TI"}} }over{M rSub{size8{0}} } 
`=` 1`-` left( 1+ { M rSub{size8{1}} }over{ M rSub{size8{0}} } right)
E rSub{size8{ital"TI"}} bold{ `=` underline{ 
	left( 1+ital"E" rSub{size8{ital"TR"}} right) `-` 
	left("1+1"right) E rSub{size8{ital"TI"}} } }
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}}`=`0 ~drarrow~ E rSub{size8{ital"TIn"}} 
`=` e rSup{size8{ -`{ ital"TIn" }over{ ital"T1" } }} 
`=` { 1+ital"E" rSub{size8{ital"TR"}} }over{ 2 }
 {} # {} # drarrow~
bold{ ital"TI" rSub{size8{n}} `=` underline{ ital"T1" left[ `"ln" left( 2 right)
	`-`"ln" left( 1+ital"E" rSub{size8{ital"TR"}} right)`right] `;~ 
   E rSub{size8{ital"TR"}} `=` e rSup{size8{ -{ ital"TReff" }over{ ital"T1" } }} } }
 {} # {} # {} # {} # alignr size8{"Q.E.D."}
} } }

........................................................................................
Slide 14: IR-Calc (2nd: Inefficient Inversion)
........................................................................................
 size 12{ alignc{ stack{
E rSub{size8{t}} `≝` e rSup{size8{-`{ t over ital"T1" }}} 
~;~ α `≝` -"cos" left( ital"FA" rSub{size8{ital"inv"}} right)
 {} # {} # 
{ M rSub{size8{1}} }over{ M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left(ital"TR"- ital"TI" right) }} 
`=` 1`-`{ E rSub{size8{ital"TR"}} }over{ E rSub{size8{ital"TI"}} } 
~;~ left( ital"TR" rightarrow ital"TR" rSub{size8{ital"eff"}} right)
 {} # {} # 
bold{S rSub{size8{0}}} `≝` { M rSub{size8{ital"TI"}} }over{M rSub{size8{0}} } 
`=` 1`-` left( 1+ital"α" { M rSub{size8{1}} }over{ M rSub{size8{0}} } right)
E rSub{size8{ital"TI"}} bold{ `=` underline{ 
	left( 1+ital"αE" rSub{size8{ital"TR"}} right) `-`
	left( 1+ital"α" right) E rSub{size8{ital"TI"}} }}
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}}`=`0 ~drarrow~ E rSub{size8{ital"TIn"}} 
`=` e rSup{size8{ -`{ ital"TIn" }over{ ital"T1" } }} 
`=` { 1+ital"αE" rSub{size8{ital"TR"}} }over{ 1+ital"α" }
 {} # {} # drarrow~
bold{ ital"TI" rSub{size8{n}} `=` underline{ ital"T1" left[ `"ln" left( 
	1+ital"α" right)`-`"ln" left( 1+ital"αE" rSub{size8{ital"TR"}} right)` right] } }
 {} # {} # {} # {} # alignr size8{"Q.E.D."}
} } }

........................................................................................
Slide 15: IR-Calc (3rd: Be T2 Prepared!)
........................................................................................
 size 12{ alignc{ stack{
E rSub{size8{t}} `≝` e rSup{size8{-`{ t over ital"T1" }}} 
~;~ α `≝` -"cos" left( ital"FA" rSub{size8{ital"inv"}} right) 
~;~ E rSub{size8{ital"T2P"}} `≝` e rSup{size8{-`{ital"T2P"}over{ital"T2"} }}
 {} # {} # 
{ M rSub{size8{1}} }over{ M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left(ital"TR"- ital"TI" right) }} 
`=` 1`-`{ E rSub{size8{ital"TR"}} }over{ E rSub{size8{ital"TI"}} } 
~;~ left( ital"TR" rightarrow ital"TR" rSub{size8{ital"eff"}} right)
 {} # {} # 
bold{S rSub{size8{0}}} `≝` { M rSub{size8{ital"TI"}} }over{M rSub{size8{0}} } 
`=` 1`-` left( 1+ital"α" { M rSub{size8{1}} }over{ M rSub{size8{0}} } right)
E rSub{size8{ital"TI"}} bold{ `=` underline{ 
	left( 1+ital"αE" rSub{size8{ital"TR"}} right)`-`
	left( 1+ital"α" right) E rSub{size8{ital"TI"}} } }
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}}`=`0 ~drarrow~ E rSub{size8{ital"TIn"}} 
`=` e rSup{size8{ -`{ ital"TIn" }over{ ital"T1" } }} 
`=` { 1+ital"αE" rSub{size8{ital"TR"}} }over{ 1+ital"α" }
 {} # {} # drarrow~
ital"TI" rSub{size8{n}} `=` underline{ ital"T1" left[ `"ln" left( 
	1+ital"α" right)`-`"ln" left( 1+ital"αE" rSub{size8{ital"TR"}} right)` right] }
 {} # {} # downarrow {} # {} # 
α rightarrow ital"αE" rSub{size8{ital"T2P"}} ~drarrow~ 
bold{ ital"TI" rSub{size8{n}} `=` underline{ ital"T1" left[ `"ln" left( 1+ital"αE"
	rSub{size8{ital"T2P"}} right)`-`"ln" left( 1+ital"αE" rSub{size8{ital"T2P"}} 
	E rSub{size8{ital"TR"}} right)` right] } }
 {} # {} # alignr size8{ital"Q.E.D."}
} } }

........................................................................................
Slide 24: DIR-Calc (Inversion juggling)
........................................................................................
 size 12{ alignc{ stack{
{ M rSub{size8{1}} }over{ M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left( ital"TR"-ital"TI" rSub{1} right) }} 
`=` 1`-`{ E rSub{size8{ital"TR"}} }over{ E rSub{size8{ital"TI" rSub{1} }} } 
~;~ left( ital"TR" rightarrow ital"TR" rSub{size8{ital"eff"}} right)
 {} # {} # 
{ M rSub{size8{2}} }over{ M rSub{size8{0}} } 
`=` 1`-`left( 1`+`α { M rSub{size8{1}} }over{ M rSub{size8{0}} } right)
E rSub{size8{ left( ital"TI" rSub{1} - ital"TI" rSub{size 8{2} } right) }} 
`=` { 1 }over{ E rSub{size8{ital"TI" rSub{2}}} } left[`E rSub{size8{ital"TI" rSub{2}}}
	+ital"αE" rSub{size8{ital"TR"}}`-`left( 1+ital"α" right)
	E rSub{size8{ital"TI" rSub{1}}}`right]
 {} # {} # 
S rSub{size8{0}} `≝` { M rSub{size8{ital"TI"}} }over{ M rSub{size8{0}} } `=` 1`-`
left( 1+ital"α" { M rSub{size8{2}} }over{ M rSub{size8{0}} } right)
	E rSub{size8{ital"TI" rSub{2}} }
 {} # {} # 
bold{ S rSub{size8{0}} `=` underline{ 1`-`left( 
	1+ital"α" right) E rSub{size8{ital"TI" rSub{2}}}`+`ital"α" left( 1+ital"α" right) 
	E rSub{size8{ital"TI" rSub{1}}}`-`
	ital"α" rSup{size8{2}} ital"E" rSub{size8{ital"TR"}} } }
 {} # {} # downarrow {} # {} # 
S rSub{size8{0}} `=` 0 ~drarrow~
 {} # {} # 
bold{ ital"TI" rSub{size8{ital"1n"}} `=` underline{ ital"T1"`"ln" left(`left[`
	ital"α" left( 1+ital"α" right)`right] `/` left[`"-1"`+`left( 1+ital"α" right)
	E rSub{size8{ital"TI" rSub{2}}}`+`
	ital"α" rSup{size8{2}} ital"E" rSub{size8{ital"TR"}}`right]`right) } }
 {} # {} # {} # alignr size8{ital"Q.E.D."}
} } }

........................................................................................
Slide 25: DIR-Calc (Oh My T2P... version)
........................................................................................
 size 12{ alignc{ stack{
{ M rSub{size8{1}} }over{ M rSub{size8{0}} } 
`=` 1`-`E rSub{size8{ left( ital"TR"-ital"TI" rSub{1} right) }} 
`=` 1`-`{ E rSub{size8{ital"TR"}} }over{ E rSub{size8{ital"TI" rSub{1} }} } 
~;~ left( ital"TR" rightarrow ital"TR" rSub{size8{ital"eff"}} right) 
~;~ ital "β"`"≝"`"αE" rSub{size8{ital"T2P"}}
 {} # {} # 
{ M rSub{size8{2}} }over{ M rSub{size8{0}} } 
`=` 1`-`left( 1`+`β { M rSub{size8{1}} }over{ M rSub{size8{0}} } right)
E rSub{size8{ left( ital"TI" rSub{1} - ital"TI" rSub{size 8{2} } right) }} 
`=` { 1 }over{ E rSub{size8{ital"TI" rSub{2}}} } left[`E rSub{size8{ital"TI" rSub{2}}}
	+ital"βE" rSub{size8{ital"TR"}}`-`left( 1+ital"β" right)
	E rSub{size8{ital"TI" rSub{1}}}`right]
 {} # {} # 
S rSub{size8{0}} `≝` { M rSub{size8{ital"TI"}} }over{ M rSub{size8{0}} } `=` 1`-`
left( 1+ital"α" { M rSub{size8{2}} }over{ M rSub{size8{0}} } right)
	E rSub{size8{ital"TI" rSub{2}} }
 {} # {} # 
bold{ S rSub{size8{0}} `=` underline{ 1`-`left( 
	1+ital"α" right) E rSub{size8{ital"TI" rSub{2}}}`+`ital"α" left( 1+ital"β" right) 
	E rSub{size8{ital"TI" rSub{1}}}`-`ital"αβE" rSub{size8{ital"TR"}} } }
 {} # {} # downarrow {} # {} # 
S rSub { size 8{0} } `=`0~ drarrow ~
 {} # {} # 
bold{ ital"TI" rSub{size8{ital"1n"}} `=` underline{ ital"T1"`"ln" left(`left[`
	ital"α" left( 1+ital"β" right)`right] `/` left[`"-1"`+`left( 1+ital"α" right)
	E rSub{size8{ital"TI" rSub{2}}}`+`ital"αβE" rSub{size8{ital"TR"}}`right]`right) } }
 {} # {} # {} # alignr size8{ital"Q.E.D."}
} } }

........................................................................................
```

