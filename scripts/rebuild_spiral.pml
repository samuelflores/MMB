
cmd.load(filename=sys.argv[1], object='spiral', format='pdb')
cmd.load(filename=sys.argv[2], object='pair', format='pdb')

resnumsA=[] ; iterate spiral and chain A and n. C5', resnumsA.append(resi)
resnumsB=[] ; iterate spiral and chain B and n. C5', resnumsB.append(resi)
resnumsB.reverse()

numres=len(resnumsA)

python
for resnr,orignr in enumerate(resnumsA):
  cmd.copy('pair'+str(resnr), 'pair')
  cmd.select('reference', '(spiral and chain A and resi '+orignr+' and (sidechain and not hydrogen)) or (spiral and chain B and resi '+str(resnumsB[resnr])+' and (sidechain and not hydrogen))')
  cmd.align('pair'+str(resnr), '(reference)')
  cmd.alter('pair'+str(resnr)+' and chain A', 'resi='+str(resnr+1))
  cmd.alter('pair'+str(resnr)+' and chain B', 'resi='+str(numres-resnr))

python end

cmd.delete('pair')
cmd.create('merged','pair*')
cmd.delete('pair*')
cmd.save('updated_'+sys.argv[1],'merged')
show sticks
zoom 'merged'
