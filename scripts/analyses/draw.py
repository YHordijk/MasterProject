import reaction


def draw(template, substituents, functional='BLYP-D3(BJ)', basis='TZ2P', numerical_quality='Good', phase='vacuum', simple=False):
	r = reaction.Reaction(template, substituents, functional=functional, basis=basis, phase=phase, numerical_quality=numerical_quality)
	print(r, r.complete, r.incomplete_reason)
	r.show(simple=simple)