
// suppress warnings in Visual studio
// TODO some of these might be legitimate
#if defined(_MSC_VER)
#pragma warning(disable:4101 4996 4018 4305 4244)
#endif

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol  ****                         *
 *                                                            *
 * molecule model module.                                     *
 *------------------------------------------------------------*/

#include "mol.h"

/*============================================================*
 *                          a t o m                           *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *              ****  mol_AtomTypeConv  ****                  *
 *                                                            *
 * convert a str atom name into an enum.                      *
 *------------------------------------------------------------*/

void 
mol_AtomTypeConv (char *s, MolAtomElementType *type)
  {

 /**************
  ***  body  ***
  **************/

  while (*s == ' ') s++;

  if (*s == 'C') {
    *type = MOL_ATOM_ELEMENT_CARBON;
    }
  else if (*s == 'N') {
    *type = MOL_ATOM_ELEMENT_NITROGEN;
    }
  else if (*s == 'O') {
    *type = MOL_ATOM_ELEMENT_OXYGEN;
    }
  else if (*s == 'H') {
    *type = MOL_ATOM_ELEMENT_HYDROGEN;
    }
  else if (*s == 'S') {
    *type = MOL_ATOM_ELEMENT_SULFUR;
    }
  else if (*s == 'P') {
    *type = MOL_ATOM_ELEMENT_PHOSPHORUS;
    }
  else {
    *type = MOL_ATOM_ELEMENT_UNKNOWN;
    }
  }

/*============================================================*
 *                m o l e c u l a r   c h a i n               *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *                 ****  mol_ChainAtomsGet  ****              *
 *                                                            *
 * get the atoms for a chain.                                 *
 *------------------------------------------------------------*/

void
mol_ChainAtomsGet (MolChain *chain, int *p_num, int **p_list) 
  {

  int num;

  int *list;

  MolChainTerm *terms;

  int num_atoms;

  MolAtom *atoms;

  char chain_id;

  int i;

 /**************
  ***  body  ***
  **************/

  chain_id = chain->id;
  mol_StructureAtomsGet (chain->struc, &num_atoms, &atoms);
  num = 0;

  for (i = 0; i < num_atoms; i++) {
    if (atoms[i].chain_id == chain_id) {
      num += 1;
      } 
    }

  if (num == 0) {
    return;
    }

  mem_Alloc (list, num, chain->model, int*); 
  num = 0;

  for (i = 0; i < num_atoms; i++) {
    if (atoms[i].chain_id == chain_id) {
      list[num] = i;
      num += 1;
      }
    }

  *p_num = num;
  *p_list = list;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_ChainCreate  ****                 *
 *                                                            *
 * create a chain.                                            *
 *------------------------------------------------------------*/

void
mol_ChainCreate (char chain_id, int num_res, int begin, int end, int rindex, 
                 MolStructure *struc, int type, MolChain **p_chain)
  {

  MolChain *chain;

  MolModel *model;

  char name[80];

 /**************
  ***  body  ***
  **************/

  mem_Alloc (chain, 1, struc->model, MolChain*);
  mol_StructureModelGet (struc, &model);

  if (type ==1 ) {
    sprintf (name, "%s:%s:%c", model->name, struc->name, chain_id);
    }
  else {
    sprintf (name, "%s:%s:solv:%c", model->name, struc->name, chain_id);
    }

  chain ->name = (char*)strdup(name);
  chain->id = chain_id;
  chain->residues.num = num_res;
  chain->residues.begin = begin;
  chain->residues.end = end;
  chain->residues.index = rindex;
  chain->residues.names = NULL;
  chain->residues.list = NULL;
  chain->type = type;

  chain->props.volume = 0.0;
  chain->props.mass = 0.0;
  chain->props.area = 0.0;

  chain->extent_set = MOL_FALSE;
  mol_ExtentSet (chain->extent, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  chain->model = NULL;
  chain->struc = struc;
  *p_chain = chain;
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  mol_ChainResiduesGet  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_ChainResiduesGet (MolChain *chain, int *p_num, MolResidue ***p_list) 
  {

  int num;

  int num_res;

  MolResidue *res, **list;

  char chain_id;

  int i, n;

  int begin, end, rindex; 

  MolResName *names;

 /**************
  ***  body  ***
  **************/

  if (chain->residues.list) { 
    *p_num = chain->residues.num;
    *p_list = chain->residues.list;
    }

  num = chain->residues.num;
  begin = chain->residues.begin;
  end = chain->residues.end;
  rindex = chain->residues.index;
  mem_Alloc (list, num, chain->model, MolResidue**);

  mol_StructureResiduesGet (chain->struc, chain->type, &num_res, &res);

  for (n = 0, i = rindex; i < rindex + num; i++, n++) {
    list[n] = res+i;
    }

  chain->residues.list = list;
  *p_num = chain->residues.num;
  *p_list = list;
  }


/*============================================================*
 *                m o l e c u l a r   m o d e l               *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *             ****  mol_MolModelCreate  ****                 *
 *                                                            *
 * create a molecule model.                                   *
 *------------------------------------------------------------*/

void
mol_MolModelCreate (const char *name, MolModel **p_model)
  {

  MolStructure *struc;

  char struct_name[500];

 /**************
  ***  body  ***
  **************/

  mem_Alloc (*p_model, 1, NULL, MolModel*);
  (*p_model)->name = (char*)strdup (name);
  (*p_model)->memory = new std::set<void*>();

  sprintf (struct_name, "%s_struct", name);
  mol_StructureCreate (struct_name, &struc, *p_model);
  struc->model = *p_model;
  mol_MolModelStructureAdd (*p_model, struc);
  (*p_model)->curr_structure = struc;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  pm_MolModelCurrStrucGet  ****         *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_MolModelCurrStrucGet (MolModel *model, MolStructure **struc) {
  *struc = model->curr_structure;
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  mol_MolModelHeaderSet  ****            *
 *                                                            *
 * add a molecule to a model.                                 *
 *------------------------------------------------------------*/

void
mol_MolModelHeaderSet (MolModel *model, char *hdr) {
  model->header = (char*)strdup(hdr);
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****  mol_MolModelStructureAdd  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_MolModelStructureAdd (MolModel *model, MolStructure *struc) {
  struc->next = model->structures;
  model->structures = struc; 
  model->num_structures += 1; 
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_MolModelStructuresGet  ****       *
 *                                                            *
 *------------------------------------------------------------*/

void 
mol_MolModelStructuresGet (MolModel *model, int *num, MolStructure **strucs) {
  *num = model->num_structures;
  *strucs = model->structures;
  }


/*============================================================*
 *                    r e s i d u e                           *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *              ****  mol_ResTypeConv  ****                   *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_ResTypeConv (char *s, MolResidueType *type)
  {

  int i;

  char tstr[4];

 /**************
  ***  body  ***
  **************/
  
  *type = MOL_RESIDUE_UNKNOWN;

  while (*s == ' ') {
    s++;
    }

  if ((s[1] == ' ') &&
      (s[2] == ' '))  { 
      s[1] = NULL;    // dropping a null here terminates the string.  This gets rid of trailing whitespace.
  }

  if (strlen(s) == 1) {
    for (i = 0; mol_res_names[i][0]; i++) {
      if (*s == *mol_res_names[i][2]) {
        *type = MolResidueType(i);
        return;
        }
      }
    }

  for (i = 0; mol_res_names[i][0]; i++) {
    if (!STRCASECMP(s, mol_res_names[i][0]) ||
        !STRCASECMP(s, mol_res_names[i][1]) || 
        !STRCASECMP(s, mol_res_names[i][2])) { 
      *type = MolResidueType(i);
      return;
      }
    }


  /*  check for gromacs modified terminal residue names.  */

  if (strlen(s) == 3) {
    tstr[0] = s[0];
    tstr[1] = s[2];
    tstr[2] = 'c';
    tstr[3] = '\0';

    for (i = 0; mol_res_names[i][0]; i++) {
      tstr[0] = mol_res_names[i][1][0];
      tstr[1] = mol_res_names[i][1][2];
      tstr[2] = 'c';
      tstr[3] = '\0';

      if (!STRCASECMP(s, tstr)) { 
        *type = MolResidueType(i);
        return;
        }

      tstr[2] = 'n';

      if (!STRCASECMP(s, tstr)) { 
        *type = MolResidueType(i);
        return;
        }
      }
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *              ****  mol_ResPropGet  ****                    *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_ResPropGet (MolResidueType type, MolResiduePropType *prop) { 
  *prop = mol_res_props[type];
  }


/*============================================================*
 *                       s t r u c t u r e                    *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureAtomAdd  ****            *
 *                                                            *
 * add an atom to a structure.                                *
 *------------------------------------------------------------*/

void
mol_StructureAtomAdd (MolStructure *struc, int het, MolAtom *atom)
  {

  MolAtomicStrucure *astruc;

  int size;

  int num;

  float xmin, xmax, ymin, ymax, zmin, zmax;

  float x, y, z;

  MolAtom *atoms;

  static float scale = 0.1f;             /*  1 angstromg = 0.1 nanometer  */

 /**************
  ***  body  ***
  **************/

  if (het) {
    astruc = &struc->het_atomic_struct;
    }
  else if (atom->res_type == MOL_RESIDUE_SOLV) {
    astruc = &struc->solvent.atomic_struct;
    }
  else {
    astruc = &struc->atomic_struct;
    }

  size = astruc->size;
  num = astruc->num_atoms;
  atoms = astruc->atoms;

  if (num == (size - 1)) {
    size += 1000;
    mem_Realloc (atoms, size, struc->model, MolAtom*);
    astruc->size = size;
    astruc->atoms = atoms;
    }


  /* note: scale atomic pos to nm.  */
  //scale = 1.0;
  x = scale*atom->pos[0];
  y = scale*atom->pos[1];
  z = scale*atom->pos[2];

  atoms[num].id = atom->id;
  atoms[num].orig_id = atom->orig_id;
  atoms[num].num = num; 
  atoms[num].name = atom->name;
  atoms[num].res_seq = atom->res_seq;
  atoms[num].insertion_code = atom->insertion_code;
  atoms[num].type = atom->type;
  atoms[num].res_type = atom->res_type;
  atoms[num].pos[0] = x;
  atoms[num].pos[1] = y;
  atoms[num].pos[2] = z; 
  atoms[num].temp = atom->temp;
  atoms[num].het = atom->het;
  atoms[num].chain_id = atom->chain_id;
  num += 1;
  astruc->num_atoms = num;

  if (het) {
    return;
    }

  if (num == 1) {
    xmin = xmax = x;
    ymin = ymax = y;
    zmin = zmax = z;
    }
  else {
    mol_ExtentGet(struc->extent, xmin, xmax, ymin, ymax, zmin, zmax);
    xmin = (x < xmin ? x : xmin);
    ymin = (y < ymin ? y : ymin);
    zmin = (z < zmin ? z : zmin);
    xmax = (x > xmax ? x : xmax);
    ymax = (y > ymax ? y : ymax);
    zmax = (z > zmax ? z : zmax);
    }

  mol_ExtentSet(struc->extent, xmin, xmax, ymin, ymax, zmin, zmax);
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****  mol_StructureAtomsGet  ****              *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureAtomsGet (MolStructure *struc, int *p_num, MolAtom **p_atoms) {
  *p_num = struc->atomic_struct.num_atoms;
  *p_atoms = struc->atomic_struct.atoms;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  pm_StructureChainAdd  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureChainAdd (MolStructure *struc, MolChain *chain, int type)
  {

 /**************
  ***  body  ***
  **************/

  if (type == 1) {
    chain->next = struc->chains;
    struc->chains = chain;
    struc->num_chains++;
    }

  else if (type == 2) {
    chain->next = struc->solvent.chains; 
    struc->solvent.chains = chain; 
    struc->solvent.num_chains++;
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****  mol_StructureChainsGet  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureChainsGet (MolStructure *struc,  int *num, MolChain **chains) {
  *num = struc->num_chains;
  *chains = struc->chains;
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****  mol_StructureChainsBuild  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureChainsBuild (MolStructure *struc, int type) 
  {

  int i, j;

  int num;

  int num_chains;

  char chain_id, chain_ids[100];

  MolChain *chains;

  int ibegin, iend, chain_index[100][4];

  int num_atoms;

  MolAtom *atoms;

  int num_res, num_chain_res, num_res_atoms;

  MolResidue *res;

  int res_seq, res_index;
 
  char insertion_code;

  MolResidueType aa_type;

  int atom_list[10000];

  MolChain *chain;

  int new_res;

  static const char *fn = "mol_StructureChainsBuild";

 /**************
  ***  body  ***
  **************/

  if (type == 1) {
    num_atoms = struc->atomic_struct.num_atoms;
    atoms = struc->atomic_struct.atoms;

/*
    if (pm_IsBase(atoms[0].aa_type)) { 
      struc->mol_type = MOL_MOLECULE_NUCLEIC_ACID; 
      }
    else {
      struc->mol_type = MOL_MOLECULE_PROTEIN; 
      }
*/
    }

  else if (type == 2) {
    num_atoms = struc->solvent.atomic_struct.num_atoms;
    atoms = struc->solvent.atomic_struct.atoms;
    }

  if (num_atoms == 0) { 
    return; 
    }

  for (i = 0; i < num_atoms; i++) {
    if (atoms[i].chain_id == ' ') { 
      atoms[i].chain_id = '0' + 1;
      }
    }

  ibegin = 0;
  chain_id = atoms[0].chain_id;
  res_seq = atoms[0].res_seq;
  insertion_code = atoms[0].insertion_code;
  num_chains = 0;
  num_res = 1;
  num_chain_res = 1;
  res_index = 0;
  
  for (i = 1; i < num_atoms; i++) {
    if (atoms[i].chain_id != chain_id) { 
      chain_ids[num_chains] = chain_id;
      chain_index[num_chains][0] = ibegin;
      chain_index[num_chains][1] = i;
      chain_index[num_chains][2] = num_chain_res;
      chain_index[num_chains][3] = res_index;

      num_chains += 1;
      num_chain_res = 1;
      ibegin = i;
      chain_id = atoms[i].chain_id;
      res_seq = atoms[i].res_seq;
      insertion_code = atoms[i].insertion_code;
      res_index = num_res;
      num_res += 1;
      }

    else if ((atoms[i].res_seq != res_seq) ||
             (atoms[i].insertion_code != insertion_code)) {
      res_seq = atoms[i].res_seq;
      insertion_code = atoms[i].insertion_code;
      num_res += 1;
      num_chain_res += 1;
      }
    }

  chain_ids[num_chains] = chain_id;
  chain_index[num_chains][0] = ibegin;
  chain_index[num_chains][1] = num_atoms-1;
  chain_index[num_chains][2] = num_chain_res;
  chain_index[num_chains][3] = res_index;
  num_chains += 1;


  /*  create residues  */

  /* note [23nov2003] change this for pdb with c-alpha atoms only.
     these guys will have only one residue per atom.                 */

  mem_Alloc (res, num_res, struc->model, MolResidue*); 
  res_seq = atoms[0].res_seq;
  insertion_code = atoms[0].insertion_code;
  aa_type = atoms[0].res_type;
  num_res = 0;
  num_res_atoms = 1;
  atom_list[0] = 0;

  for (i = 1; i < num_atoms; i++) {
    if ((atoms[i].res_seq != res_seq) ||
        (atoms[i].insertion_code != insertion_code)) { 

      res[num_res].id = res_seq; 
      res[num_res].insertion_code = insertion_code;
      res[num_res].atoms.num = num_res_atoms; 
      res[num_res].type = aa_type; 
      mem_Alloc (res[num_res].atoms.list, num_res_atoms, struc->model, int*); 
      res[num_res].active = MOL_TRUE; 

      for (j = 0; j < num_res_atoms; j++) {
        res[num_res].atoms.list[j] = atom_list[j]; 
        }

      res_seq = atoms[i].res_seq;
      insertion_code = atoms[i].insertion_code;
      aa_type = atoms[i].res_type;
      num_res += 1;
      num_res_atoms = 0;
      }

    atom_list[num_res_atoms] = i;
    num_res_atoms += 1;
    }

  res[num_res].id = res_seq; 
  res[num_res].insertion_code = insertion_code;
  res[num_res].atoms.num = num_res_atoms; 
  res[num_res].type = aa_type; 
  mem_Alloc (res[num_res].atoms.list, num_res_atoms, struc->model, int*); 

  for (j = 0; j < num_res_atoms; j++) {
    res[num_res].atoms.list[j] = atom_list[j]; 
    }

  num_res += 1;


  /*  create chains  */

  for (i = 0; i < num_chains; i++) {
    chain_id = chain_ids[i];

    if (chain_id == ' ') {
      chain_id = '0' + 1; 
      chain_ids[i] = chain_id;
      }

    ibegin = chain_index[i][0];
    iend = chain_index[i][1];
    res_index = chain_index[i][3];
    num_chain_res = chain_index[i][2];
    mol_ChainCreate (chain_id, num_chain_res, ibegin, iend, res_index, struc, type, 
                    &chain);
    mol_StructureChainAdd (struc, chain, type);
    }

  if (type == 1) {
    struc->num_res = num_res;
    struc->residues = res;
    }

  else if (type == 2) {
    struc->solvent.num_res = num_res;
    struc->solvent.residues = res;
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  mol_StructureCreate  ****              *
 *                                                            *
 * create a structure object.                                 *
 *------------------------------------------------------------*/

void
mol_StructureCreate (char *name, MolStructure **p_struc, MolModel *model)
  {

  int i, size;

  MolStructure *struc;

 /**************
  ***  body  ***
  **************/

  mem_Alloc (struc, 1, model, MolStructure*);
  struc ->name = (char*)strdup(name);

  size = 1000;
  struc->atomic_struct.size = size;
  mem_Alloc (struc->atomic_struct.atoms, size, model, MolAtom*);

  size = 1000;
  struc->het_atomic_struct.size = size;
  mem_Alloc (struc->het_atomic_struct.atoms, size, model, MolAtom*);

  size = 1000;
  struc->solvent.atomic_struct.size = size;
  mem_Alloc (struc->solvent.atomic_struct.atoms, size, model, MolAtom*);

  for (i = 0; i < MOL_ATOM_ELEMENT_MAX; i++) {
    struc->atomic_struct.atom_groups[i] = NULL; 
    struc->het_atomic_struct.atom_groups[i] = NULL;
    struc->solvent.atomic_struct.atom_groups[i] = NULL; 
    }

  size = 10;
  struc->secondary.helix.size = size;
  mem_Alloc (struc->secondary.helix.list, size, model, MolHelix*); 

  size = 10;
  struc->secondary.sheet.size = size;
  mem_Alloc (struc->secondary.sheet.list, size, model, MolSheet*); 

  *p_struc = struc;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureHelixAdd  ****           *
 *                                                            *
 * add a helix to a structure.                                *
 *------------------------------------------------------------*/

void
mol_StructureHelixAdd (MolStructure *struc, MolHelix *helix)
  {

  MolStructureSecHelix *hstruc;

  MolHelix *hlist;

  int size;

  int num;

  float x, y, z;

  MolAtom *atoms;

 /**************
  ***  body  ***
  **************/

  hstruc = &struc->secondary.helix;
  num = hstruc->num;
  size = hstruc->size;
  hlist = hstruc->list;

  if (num == (size - 1)) {
    size += 1000;
    mem_Realloc (hlist, size, struc->model, MolHelix*);
    hstruc->size = size;
    hstruc->list = hlist;
    }

  memcpy (hlist+num, helix, sizeof(MolHelix));
  num += 1;
  hstruc->num = num;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureHelicesGet  ****         *
 *                                                            *
 * get helices for a structure.                               *
 *------------------------------------------------------------*/

void
mol_StructureHelicesGet (MolStructure *struc, int *num, MolHelix **list)
  {

  MolStructureSecHelix *hstruc;

 /**************
  ***  body  ***
  **************/

  hstruc = &struc->secondary.helix;
  *num = hstruc->num;
  *list = hstruc->list;
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****  pm_StructureResiduesGet  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void 
mol_StructureResiduesGet (MolStructure *struc, int type, int *num_res, MolResidue **res)
  {

  char id;

 /**************
  ***  body  ***
  **************/

  *num_res = 0;

  if (type == 1) {
    *num_res = struc->num_res;
    *res = struc->residues;
    }

  else if (type == 2) {
    *num_res = struc->solvent.num_res;
    *res = struc->solvent.residues;
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureSheetAdd  ****           *
 *                                                            *
 * add a sheet to a structure.                                *
 *------------------------------------------------------------*/

void
mol_StructureSheetAdd (MolStructure *struc, MolSheet *sheet)
  {

  MolStructureSecSheet *sstruc;

  MolSheet *list;

  int size;

  int num;

 /**************
  ***  body  ***
  **************/

  sstruc = &struc->secondary.sheet;
  num = sstruc->num;
  size = sstruc->size;
  list = sstruc->list;

  if (num == (size - 1)) {
    size += 1000;
    mem_Realloc (list, size, struc->model, MolSheet*);
    sstruc->size = size;
    sstruc->list = list;
    }

  memcpy (list+num, sheet, sizeof(MolSheet));
  num += 1;
  sstruc->num = num;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureSheetsGet  ****          *
 *                                                            *
 * get sheets for a structure.                                *
 *------------------------------------------------------------*/

void
mol_StructureSheetsGet (MolStructure *struc, int *num, MolSheet **list)
  {

  MolStructureSecSheet *sstruc;

 /**************
  ***  body  ***
  **************/

  sstruc = &struc->secondary.sheet;
  *num = sstruc->num;
  *list = sstruc->list;
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  mol_StructureTermAdd  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureTermAdd (MolStructure *struc, MolChainTerm *term, MolModel *model)
  {

  MolChainTerm *tp;

 /**************
  ***  body  ***
  **************/

  mem_Alloc (tp, 1, model, MolChainTerm*); 
  tp->id = term->id;
  tp->chain_id = term->chain_id;
  tp->res_seq = term->res_seq;
  tp->insertion_code = term->insertion_code;
  strcpy (tp->res_name, term->res_name);
  tp->next = struc->terms;
  struc->terms = tp;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_StructureModelGet  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_StructureModelGet (MolStructure *struc, MolModel **model) { 
  *model = struc->model; 
  }


/*============================================================*
 *                       d a t a b a s e                      *
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_DbFormatConv  ****                *
 *                                                            *
 *------------------------------------------------------------*/

void 
mol_DbFormatConv (char *format,  MolDbType *type)
  {

 /**************
  ***  body  ***
  **************/

  if (!strcmp(format, "pdb")) {
    *type = MOL_DB_PDB;
    }

  else if (!strcmp(format, "gromacs")) {
    *type = MOL_DB_GROMACS;
    }

  else {
   *type = MOL_DB_UNKNOWN;
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol_DbRead  ****                  *
 *                                                            *
 * read a molecule database.                                  *
 *------------------------------------------------------------*/

void 
mol_DbRead (const char *name, const char *file_name, MolDbType db_type, MolModel **p_model) 
  {

  FILE *fp;

  char line[1000];

  int n, i;

  void (*func)(FILE *fp, char *line, MolModel *model);

  MolModel *model;

  MolStructure *struc;

 /**************
  ***  body  ***
  **************/

  fp = fopen (file_name, "r");

  if (!fp) {
    fprintf (stderr, "\n  ****  error: can't open file [%s] \n", file_name);
    return;
    }

  mol_MolModelCreate (name, &model);

  if (db_type == MOL_DB_GROMACS) {
    mol_DbGromacsRead (model, fp);
    }

  else if (db_type == MOL_DB_PDB) {
    while (1) {
      if (!mol_DbLineGet(fp, line)) {
        break;
        }

      for (i = 0; pdb_table[i].name != NULL; i++) {
        n = strlen (pdb_table[i].name);
  
        if (!STRNCASECMP(line, pdb_table[i].name, n)) {
          func = pdb_table[i].func;
          (*func)(fp, line, model);
          }
        }
      }
    }

  mol_MolModelCurrStrucGet (model, &struc);

  /* build protein chains */

  mol_StructureChainsBuild (struc, 1); 


  /* build solvent chains */

  mol_StructureChainsBuild (struc, 2); 

  *p_model = model;
  fclose (fp);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol_DbGromacsRead  ****           *
 *                                                            *
 * read a gromacs molecule database.                          *
 *------------------------------------------------------------*/

void 
mol_DbGromacsRead (MolModel *model, FILE *fp) 
  {

  char line[1000];

  int i, j, n;

  MolStructure *struc;

  int num_atoms;

  int atom_id, res_id;

  char *s, c, res_name[10], atom_name[10];

  float x, y, z, vx, vy, vz;

  MolAtom atom;

  static const char *fmt = "%5d %5s %5s %5d %8.3f% 8.3f %8.3f ";

 /**************
  ***  body  ***
  **************/

  char* retp; // to avoid warnings
  retp=fgets (line, 1000, fp);
  retp=fgets (line, 1000, fp);
  num_atoms = atoi (line);
  mol_MolModelCurrStrucGet (model, &struc);

  for (i = 0; i < num_atoms; i++) {
    retp = fgets (line, 1000, fp);
    mol_DbRecordIntGet (line, 1, 5, &res_id);
    mol_DbRecordStrGet (line, 6, 10, res_name);
    mol_DbRecordStrGet (line, 11, 15, atom_name);
    mol_DbRecordIntGet (line, 16, 20, &atom_id);
    mol_DbRecordRealGet (line, 21, 28, &x);
    mol_DbRecordRealGet (line, 29, 36, &y);
    mol_DbRecordRealGet (line, 37, 44, &z);


    /*  modify atom name to conform with pdb */

    s = atom_name; 
    while (*s == ' ') s++;
    sprintf (atom_name, " %s", s); 

    for (j = strlen(atom_name); j < 4; j++) {
      atom_name[j] = ' ';
      }

    atom_name[4] = '\0';


    /* coords in nm so scale to angstroms.  braindead.  */

    x *= 10.0;
    y *= 10.0;
    z *= 10.0;


    /* get atom element  */

    for (j = 0; j < strlen(atom_name); j++) {
      if (atom_name[j] != ' ') {
        mol_AtomTypeConv (atom_name+j, &atom.type);
        break;
        }
      }

    atom.name = (char *)strdup(atom_name);
    atom.id = atom_id; 
    atom.orig_id = atom_id; 
    atom.res_seq = res_id; 
    //atom.insertion_code = insertion_code; 
    atom.het = MOL_FALSE;
    atom.chain_id = 'A'; 
    mol_PointSet3 (atom.pos, x, y, z); 
    res_name[3] = '\0';
    mol_ResTypeConv (res_name, &atom.res_type);
    mol_ResPropGet (atom.res_type, &atom.res_prop);


    /*  add the atom to the model.  */

    mol_StructureAtomAdd (struc, MOL_FALSE, &atom);
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                     ****  mol_DbLineGet  ****              *
 *                                                            *
 *------------------------------------------------------------*/

int
mol_DbLineGet (FILE *fp, char *line)
  {

 /**************
  ***  body  ***
  **************/

  while (1) {
    if (fgets(line, 2000, fp) == NULL) {
      return (0);
      }
    else {
      mol_DbLineStrip (line);

      if (line[0] != '\n') {
        return (1);
        }
      }
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol_DbLineStrip  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbLineStrip (char *s)
  {

  int i, n;

 /**************
  ***  body  ***
  **************/

  for (i = 0, n = 0; i < strlen(s); i++) {
    if (n != 0) {
      s[n++] = s[i];
      }
    else if (s[i] != ' ') {
      s[n++] = s[i];
      }
    }

  s[n] = '\0';
  }


/*------------------------------------------------------------*
 *                                                            *
 *              ****  mol_DbPdbHeaderProc  ****               *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbHeaderProc (FILE *fp, char *line, MolModel *model)
  {

  char *s;

 /**************
  ***  body  ***
  **************/

  s = (char *)strtok (line, " ");
  s = (char *)strtok (NULL, " ");
  mol_MolModelHeaderSet (model, s);
  }


/*------------------------------------------------------------*
 *                                                            *
 *              ****  mol_DbPdbCompProc  ****                 *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbCompProc (FILE *fp, char *line, MolModel *model)
  {

  char *s;

 /**************
  ***  body  ***
  **************/

  //s = (char *)strtok (line, " ");
  //mol_MolModelHeaderSet (model, s);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbPdbHelixProc  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbHelixProc (FILE *fp, char *line, MolModel *model)
  {

  char str[80];

  MolStructure *struc;

  MolHelix helix;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 8, 10, &helix.num);
  mol_DbRecordStrGet (line, 12, 14, helix.id);

  mol_DbRecordStrGet (line, 16, 18, helix.init_res_name);
  mol_DbRecordStrGet (line, 20, 20, &helix.init_chain_id);
  mol_DbRecordIntGet (line, 22, 25, &helix.init_seq_num);

  mol_DbRecordStrGet (line, 28, 30, helix.term_res_name);
  mol_DbRecordStrGet (line, 32, 32, &helix.term_chain_id);
  mol_DbRecordIntGet (line, 34, 37, &helix.term_seq_num);

  mol_DbRecordIntGet (line, 39, 40, &helix.helixClass);
  mol_DbRecordIntGet (line, 72, 76, &helix.length);

  mol_MolModelCurrStrucGet (model, &struc);
  mol_StructureHelixAdd (struc, &helix);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_DbPdbSeqresProc  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbSeqresProc (FILE *fp, char *line, MolModel *model)
  {

  MolStructure *struc;

  int res_seq;
  
  char insertion_code;

  int num_res, n;

  char chain_id;

  char *s, str[80];

  MolResName *residues;

  MolChain *chain;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 9,  10, &res_seq);
  mol_DbRecordIntGet (line, 14, 17, &num_res);
  chain_id = line[11];
  mem_Alloc (residues, num_res, model, MolResName*);
  n = 0;

  while (1) {
    mol_DbRecordStrGet (line, 20, 70, str);
    s = (char *)strtok (str, " ");

    while (s && (*s != '\n') && (*s != ' ')) {
      strcpy (residues[n++], s);
      s = (char *)strtok (NULL, " ");
      }

    if (n == num_res) {
      break;
      }

    if (!mol_DbLineGet(fp, line)) {
      break;
      }
    }

  /*
  mol_MolModelCurrStrucGet (model, &struc);
  mol_ChainCreate (chain_id, num_res, residues, model, struc, &chain);
  mol_StructureChainAdd (struc, chain); 
  */
  }


/*------------------------------------------------------------*
 *                                                            *
 *                 ****  mol_DbPdbHetAtomProc  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbHetAtomProc (FILE *fp, char *line, MolModel *model)
  {

  char str[80];

  MolStructure *struc;

  MolAtom atom;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 7, 11, &atom.id);
  mol_DbRecordStrGet (line, 13, 16, str);
  atom.name = (char *)strdup(str);

  mol_DbRecordIntGet (line, 23, 26, &atom.res_seq);
  // 27 is insertion code
  atom.insertion_code = line[26]; // column 27 is the insertion code
  mol_DbRecordRealGet (line, 31, 38, &atom.pos[0]);
  mol_DbRecordRealGet (line, 39, 46, &atom.pos[1]);
  mol_DbRecordRealGet (line, 47, 54, &atom.pos[2]);

  mol_DbRecordRealGet (line, 61, 66, &atom.temp);
  mol_DbRecordStrGet (line, 77, 78, str);

  /* atom element is always the 2nd char  */
  mol_AtomTypeConv (str+2, &atom.type);
  atom.het = MOL_TRUE;
  atom.chain_id = line[21];


  /*  add the atom to the model.  */

  mol_MolModelCurrStrucGet (model, &struc);
  mol_StructureAtomAdd (struc, MOL_TRUE, &atom);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                     ****  mol_DbPdbTerProc  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbTerProc (FILE *fp, char *line, MolModel *model)
  {

  MolStructure *struc;

  MolChainTerm term;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 7, 11, &term.id);
  mol_DbRecordStrGet (line, 18, 20, term.res_name);
  term.chain_id = line[21];
  mol_DbRecordIntGet (line, 23, 26, &term.res_seq);

  term.insertion_code = line[26]; // column 27 is the insertion code
  mol_MolModelCurrStrucGet (model, &struc);
  mol_StructureTermAdd (struc, &term, model);
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  mol_DbPdbSheetProc  ****               *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbSheetProc (FILE *fp, char *line, MolModel *model)
  {

  char str[80];

  MolStructure *struc;

  MolSheet sheet;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 8, 10, &sheet.num);
  mol_DbRecordStrGet (line, 12, 14, sheet.id);

  mol_DbRecordIntGet (line, 15, 16, &sheet.num_strands);
  mol_DbRecordStrGet (line, 18, 20, sheet.init_res_name);

  mol_DbRecordStrGet (line, 22, 22, &sheet.init_chain_id);
  mol_DbRecordIntGet (line, 23, 26, &sheet.init_seq_num);
  mol_DbRecordStrGet (line, 27, 27, &sheet.init_icode);

  mol_DbRecordStrGet (line, 29, 31, sheet.term_res_name);
  mol_DbRecordStrGet (line, 33, 33, &sheet.term_chain_id);
  mol_DbRecordIntGet (line, 34, 37, &sheet.term_seq_num);
  mol_DbRecordStrGet (line, 38, 38, &sheet.end_icode);

  mol_DbRecordIntGet (line, 39, 40, &sheet.sense);


  mol_MolModelCurrStrucGet (model, &struc);
  mol_StructureSheetAdd (struc, &sheet);

  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbPdbAtomProc  ****             *
 *                                                            *
 * process pdb atom line.                                     *
 *------------------------------------------------------------*/

void
mol_DbPdbAtomProc (FILE *fp, char *line, MolModel *model)
  {

  char str[80], aa_name[4];

  int n, i;

  MolAtom atom;

  MolStructure *struc;

 /**************
  ***  body  ***
  **************/

  mol_DbRecordIntGet (line, 7, 11, &atom.id);
  atom.orig_id = atom.id;
  mol_DbRecordStrGet (line, 13, 16, str);
  atom.name = (char *)strdup(str); 
  mol_DbRecordStrGet (line, 18, 20, aa_name);
  mol_ResTypeConv (aa_name, &atom.res_type);
  mol_ResPropGet (atom.res_type, &atom.res_prop);

  /* alternate location indicator. i guess use the 1st one.  */

  if ((line[16] != 'A') && (line[16] != ' ')) {
    return;
    }

  mol_DbRecordIntGet (line, 23, 26, &atom.res_seq);
  atom.insertion_code = line[26]; // column 27 is the insertion code
  mol_DbRecordRealGet (line, 31, 38, &atom.pos[0]);
  mol_DbRecordRealGet (line, 39, 46, &atom.pos[1]);
  mol_DbRecordRealGet (line, 47, 54, &atom.pos[2]);

  mol_DbRecordRealGet (line, 61, 66, &atom.temp);
  mol_DbRecordStrGet (line, 77, 78, str);
  /*
  mol_AtomTypeConv (str, &atom.type);
  */

  if (atom.name[0] != ' ') {
    mol_AtomTypeConv (&atom.name[0], &atom.type);
    }
  else {
    mol_AtomTypeConv (&atom.name[1], &atom.type);
    }

  atom.het = MOL_FALSE;
  atom.chain_id = line[21];


  /*  add the atom to the molecule.  */

  mol_MolModelCurrStrucGet (model, &struc);
  mol_StructureAtomAdd (struc, MOL_FALSE, &atom);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbPdbConectProc  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbPdbConectProc (FILE *fp, char *line, MolModel *mol)
  {

  char *s;

  int n;

  int conn[1000];

 /**************
  ***  body  ***
  **************/

  return;
  s = strtok (line, " ");
  s = strtok (NULL, " ");
  n = 0;

  while (*s != '\n') {
    conn[n++] = atoi (s) - 1;
    s = strtok (NULL, " ");
    }

/*
  mol_MolAtomConnAdd (mol, n, conn);
*/
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbRecordIntGet  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbRecordIntGet (char *line, int start, int end, int *val)
  {

  char str[80];

 /**************
  ***  body  ***
  **************/

  mol_DbRecordStrGet (line, start, end, str);
  *val = atoi (str);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbStrParse  ****                *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbStrParse (char *line, char *val)
  {

  int i, j;

 /**************
  ***  body  ***
  **************/

  for (i = 0; i < strlen(line); ++i) {
    if (line[i] == '"') {
      break;
      }
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbRecordStrGet  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbRecordStrGet (char *line, int start, int end, char *str)
  {

  int i, j;

 /**************
  ***  body  ***
  **************/

  start--;
  end--;

  for (j = 0, i = start; i <= end; i++, j++) {
    str[j] = line[i];
    }

  if (start != end) {
    str[j] = '\0';
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_DbRecordRealGet  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void
mol_DbRecordRealGet (char *line, int start, int end, float *val)
  {

  char str[80];

 /**************
  ***  body  ***
  **************/

  mol_DbRecordStrGet (line, start, end, str);
  *val = atof (str);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol_DbMolSurfRead  ****           *
 *                                                            *
 * read a pm surface file.                                    *
 *------------------------------------------------------------*/

void 
mol_DbMolSurfRead (char *file_name, MolSurface *surf) 
  {

  FILE *fp;

  char line[1000];

  int n, i;

 /**************
  ***  body  ***
  **************/

  fp = fopen (file_name, "r");

  if (!fp) {
    fprintf (stderr, "\n  ****  error: can't open file [%s] \n", file_name);
    return;
    }

  mol_DbMolSurfParse (fp, surf);
  fclose (fp);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  mol_DbMolSurfParse  ****          *
 *                                                            *
 * read a pm surface object.                                  *
 *------------------------------------------------------------*/

void 
mol_DbMolSurfParse (FILE *fp, MolSurface *surf)
  {

  int num_verts;

  MolPoint3 *verts;

  int num_tri;

  MolIntPoint3 *conn;

  char *s, line[1000];

  char name[1000];

  int n, i, id;

 /**************
  ***  body  ***
  **************/

  while (1) {
    if (!mol_DbLineGet(fp, line)) {
      break;
      }

    if (!strncmp(line, "surface", 7)) {
      s = strtok (line, "\"");
      s = strtok (NULL, "\"");
      strcpy (name, s);
      }

    else if (!strncmp(line, "binary", 6)) {
      mol_DbMolSurfBinRead (fp, surf);
      return;
      }

    else if (!strncmp(line, "number of vertices", 18)) {
      sscanf (line, "number of vertices %d \n", &num_verts);
      mem_Alloc (verts, num_verts, NULL, MolPoint3*);
      mol_DbLineGet (fp, line);

      for (i = 0; i < num_verts; i++) {
        int nn=fscanf (fp, "%d %f %f %f \n", &id, &verts[i][0], &verts[i][1], &verts[i][2]);
        verts[i][0] *= 0.1;
        verts[i][1] *= 0.1;
        verts[i][2] *= 0.1;
        }

      mol_DbLineGet (fp, line);
      }

    else if (!strncmp(line, "number of triangles", 19)) {
      sscanf (line, "number of triangles %d \n", &num_tri);
      mem_Alloc (conn, num_tri, NULL, MolIntPoint3*);
      mol_DbLineGet (fp, line);

      for (i = 0; i < num_tri; i++) {
        int nn=fscanf (fp, "%d %d %d %d \n", &id, &conn[i][0], &conn[i][1], &conn[i][2]);
        }

      mol_DbLineGet (fp, line);
      }

    else if (!strncmp(line, "surface end", 11)) {
      break;
      }
    }

  surf->num_verts = num_verts; 
  surf->verts = verts; 
  surf->num_tri = num_tri; 
  surf->conn = conn; 
  }


/*------------------------------------------------------------*
 *                                                            *
 *             ****   mol_DbMolSurfBinRead  ****              *
 *                                                            *
 * read a pm surface object.                                  *
 *------------------------------------------------------------*/

void
mol_DbMolSurfBinRead (FILE *fp, MolSurface *surf)
  {

  int num_verts;

  MolPoint3 *verts;

  int num_tri;

  MolIntPoint3 *conn;

  int i;

 /**************
  ***  body  ***
  **************/

  size_t nRead; // to avoid warnings
  nRead=fread (&num_verts, sizeof(int), 1, fp);
  nRead=fread (&num_tri, sizeof(int), 1, fp);
  mem_Alloc (verts, num_verts, NULL, MolPoint3*);
  nRead=fread (verts, sizeof(MolPoint3), num_verts, fp);

  mem_Alloc (conn, num_tri, NULL, MolIntPoint3*);
  nRead=fread (conn, sizeof(MolIntPoint3), num_tri, fp);

  surf->num_verts = num_verts; 
  surf->verts = verts; 
  surf->num_tri = num_tri; 
  surf->conn = conn; 

  for (i = 0; i < num_verts; i++) {
    verts[i][0] *= 0.1;
    verts[i][1] *= 0.1;
    verts[i][2] *= 0.1;
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                ****  mol_MemCallocFunc  ****               *
 *                                                            *
 *------------------------------------------------------------*/

void *
mol_MemCallocFunc (const char *name, const char *file, int line, unsigned num, unsigned size, MolModel* model)
  {

  void *p;

  int pid;

 /**************
  ***  body  ***
  **************/

  p = NULL;

  if ((num == 0) || (size == 0)) {
    p = NULL;
    }

  else {
    p = (void*)calloc((size_t)num, (size_t)size);

    if (p == NULL) {
      mol_msg ("\n**** error: calloc failed. \n", 1);
      }
    }

  if (p != NULL && model != NULL)
    model->memory->insert(p);
  return p;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  mol_MemReallocFunc  ****            *
 *                                                            *
 *------------------------------------------------------------*/

void *
mol_MemReallocFunc (const char *name, const char *file, int line, void *ptr, unsigned size, MolModel* model)
  {

  void *p;

  int pid;

 /**************
  ***  body  ***
  **************/

  p = NULL;

  if (ptr == 0) {
    p = (void*)malloc((size_t)size);;
    mol_msg ("\n**** waring: realloc null ptr. \n", 1);
    }

  else {
    if (model != NULL)
      model->memory->erase(model->memory->find(ptr));
    p = (void*)realloc(ptr, (size_t)size);

    if (p == NULL) {
      mol_msg ("\n**** error: realloc failed. \n", 1);
      }
    }

  if (p != NULL && model != NULL)
    model->memory->insert(p);
  return p;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                   ****  mol_MemFreeModel  ****             *
 *                                                            *
 * clean up memory that was allocated for the model           *
 *------------------------------------------------------------*/


void mol_MemFreeModel(MolModel* model) {
    for (std::set<void*>::const_iterator iter = model->memory->begin(); iter != model->memory->end(); ++iter)
        free(*iter);
    delete model->memory;
    free(model);
}



/*------------------------------------------------------------*
 *                                                            *
 *                      ****  msg  ****                       *
 *                                                            *
 * output a message to stderr.                                *
 *------------------------------------------------------------*/

void
mol_msg (const char *format, ...)
  {

  va_list ap;

  int c;

  int ia;

  int i;

  float fa;

  char *sa, *fn;

  char des[3], msg[1000], gui_msg[1000], str[1000], tf[1000];

  int n;

  int pr;

 /**************
  ***  body  ***
  **************/

  i = 0;
  n = 0;
  va_start (ap, format);
  strcpy (tf, format);
  pr = va_arg (ap, int); 

  if (!pr) {
    return;
    }

  msg[0] = '\0';


  /*  parse the message format. when a 
      format (%) character is encountered
      extract an argument from the arglist
      of this function and print it into
      msg[] using the format parsed so far.  */
      
  while ((c = *format++) != 0) {
    if (i == 0) {
      if (c == '%') {
        des[i++] = c;
        }
      tf[n++] = c;
      }
    else if (i == 1) {
      des[i++] = c;
      des[i] = '\0';
      tf[n++] = c;
      tf[n] = '\0';


      /*  int format.  */

      if (!strcmp(des, "%d")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  float format.  */

      else if (!strcmp(des, "%f")) {
	fa = va_arg (ap, double);
        sprintf (str, tf, fa);
        sprintf (msg, "%s%s", msg, str);
	}


     /*  float format.  */

     else if (!strcmp(des, "%g")) {
       fa = va_arg (ap, double);
       sprintf (str, tf, fa);
       sprintf (msg, "%s%s", msg, str);
       }


      /*  hex format.  */

      else if (!strcmp(des, "%x")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  string format.  */

      else if (!strcmp(des, "%s")) {
	sa = va_arg (ap, char *);
        sprintf (str, tf, sa);
        sprintf (msg, "%s%s", msg, str);
	}

      i = 0;
      n = 0;
      }
    }

  va_end (ap);
  tf[n] = '\0';

  fprintf (stderr, "%s%s ", msg, tf); 
  fflush (stderr);
  }


