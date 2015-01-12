/**
 * This file and mol.cpp contain code for parsing PDB files (among others).  It was written
 * by David Parker, and modified by Peter Eastman.
 */


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  mol  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

#include <set>
#include <cstdio>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

// <strings.h>, strcasecmp, and strncasecmp are not POSIX standard
// this stanza attempts to invoke suitable replacements for windows
#ifdef _WIN32
#define STRCASECMP stricmp
#define STRNCASECMP strnicmp
#else
#include <strings.h>
#define STRCASECMP strcasecmp
#define STRNCASECMP strncasecmp
#endif



typedef  float  (MolPoint3)[3];
typedef  float  (MolVector3)[3];
typedef  float  (MolExtent)[2][3];
typedef  int    (MolIntPoint3)[3];

#define mol_PointSet3(point, v1, v2, v3) \
          point[0] = v1;                 \
          point[1] = v2;                 \
          point[2] = v3                  \

#define mol_ExtentGet(extent, xmin, xmax, ymin, ymax, zmin, zmax) \
          xmin = extent[0][0];                                    \
          xmax = extent[1][0];                                    \
          ymin = extent[0][1];                                    \
          ymax = extent[1][1];                                    \
          zmin = extent[0][2];                                    \
          zmax = extent[1][2]                                     \

#define mol_ExtentSet(extent, xmin, xmax, ymin, ymax, zmin, zmax) \
          extent[0][0] = xmin;                                    \
          extent[1][0] = xmax;                                    \
          extent[0][1] = ymin;                                    \
          extent[1][1] = ymax;                                    \
          extent[0][2] = zmin;                                    \
          extent[1][2] = zmax

#define MOL_TRUE  1
#define MOL_FALSE 0

typedef enum MolDbType {
  MOL_DB_UNKNOWN,
  MOL_DB_PDB,
  MOL_DB_GROMACS,
  } MolDbType;

typedef enum MolMoleculeType {
  MOL_MOLECULE_UNKNOWN,
  MOL_MOLECULE_RNA,
  MOL_MOLECULE_DNA,
  MOL_MOLECULE_NUCLEIC_ACID,
  MOL_MOLECULE_PROTEIN,
  } MolMoleculeType;

typedef struct MolSurface {
  int num_verts;
  MolPoint3 *verts;
  int num_tri;
  MolIntPoint3 *conn;
  } MolSurface;

typedef enum MolResiduePropType {
  MOL_RESIDUE_PROP_UNKNOWN,
  MOL_RESIDUE_PROP_ACIDIC,
  MOL_RESIDUE_PROP_BASIC,
  MOL_RESIDUE_PROP_HYDROPHOBIC,
  MOL_RESIDUE_PROP_POLAR
  } MolResiduePropType;

typedef enum MolResidueType {
  MOL_RESIDUE_UNKNOWN,
  MOL_RESIDUE_ALANINE,
  MOL_RESIDUE_ARGININE,
  MOL_RESIDUE_ASPARAGINE,
  MOL_RESIDUE_ASPARTIC_ACID,
  MOL_RESIDUE_CYSTEINE,
  MOL_RESIDUE_GLUTAMIC_ACID,
  MOL_RESIDUE_GLUTAMINE,
  MOL_RESIDUE_GLYCINE,
  MOL_RESIDUE_HISTIDINE,
  MOL_RESIDUE_ISOLEUCINE,
  MOL_RESIDUE_LEUCINE,
  MOL_RESIDUE_LYSINE,
  MOL_RESIDUE_METHIONINE,
  MOL_RESIDUE_PHENYLALANINE,
  MOL_RESIDUE_PROLINE,
  MOL_RESIDUE_SERINE,
  MOL_RESIDUE_THREONINE,
  MOL_RESIDUE_TRYPTOPHAN,
  MOL_RESIDUE_TYROSINE,
  MOL_RESIDUE_VALINE,
  MOL_RESIDUE_SOLV,
  MOL_RESIDUE_ADP,
  MOL_RESIDUE_CHLORINE,
  MOL_RESIDUE_HOH,
  MOL_RESIDUE_ADENOSINE,
  MOL_RESIDUE_GUANOSINE,
  MOL_RESIDUE_CYTOSINE,
  MOL_RESIDUE_URIDINE,
  MOL_RESIDUE_THYMINE,
  MOL_RESIDUE_URIDINE2,
  MOL_RESIDUE_CYTOSINE2,
  MOL_RESIDUE_GUANOSINE2,
  MOL_RESIDUE_URIDINE3,
  MOL_RESIDUE_CYTOSINE3,
  MOL_RESIDUE_GUANOSINE3,
  MOL_RESIDUE_URIDINE4,
  MOL_RESIDUE_ADENOSINE2,
  MOL_RESIDUE_GUANOSINE4,
  MOL_RESIDUE_GUANOSINE5,
  MOL_RESIDUE_MAX_NUM
  } MolResidueType;

static const char *mol_res_names[MOL_RESIDUE_MAX_NUM+1][3] = {
  {"unknown",       "unknown", "unknown"},
  {"alanine",       "ala",     "a"},
  {"arginine",      "arg",     "r"},
  {"asparagine",    "asn",     "n"},
  {"aspartic acid", "asp",     "d"},
  {"cysteine",      "cys",     "c"},
  {"glutamic acid", "glu",     "e"},
  {"glutamine",     "gln",     "q"},
  {"glycine",       "gly",     "g"},
  {"histidine",     "his",     "h"},
  {"isoleucine",    "ile",     "i"},
  {"leucine",       "leu",     "l"},
  {"lysine",        "lys",     "k"},
  {"methionine",    "met",     "m"},
  {"phenylalanine", "phe",     "f"},
  {"proline",       "pro",     "p"},
  {"serine",        "ser",     "s"},
  {"threonine",     "thr",     "t"},
  {"tryptophan",    "trp",     "w"},
  {"tyrosine",      "tyr",     "y"},
  {"valine",        "val",     "v"},
  {"solvent",       "sol",     "h"},
  {"adp",           "adp",     "h"},
  {"cl",            "cli",     "h"},
  {"hoh",           "hoh",     "h"},
  {"adenosine",     "adn",     "A"},
  {"guanosine",     "gua",     "G"},
  {"cytosine",      "cyt",     "C"},
  {"uridine",       "ura",     "U"},
  {"thymine",       "thy",     "T"},
  {"uridine",       "h2u",     "U"},
  {"cytosine",      "omc",     "C"},
  {"guanosine",     "omg",     "G"},
  {"uridine",       "psu",     "U"},
  {"cytosine",      "5mc",     "C"},
  {"guanosine",     "7mg",     "G"},
  {"uridine",       "5mu",     "U"},
  {"adenosine",     "1ma",     "A"},
  {"guanosine",     "2mg",     "G"},
  {"guanosine",     "m2g",     "G"},
  {"",              "",        ""}};


 static MolResiduePropType mol_res_props[MOL_RESIDUE_MAX_NUM] = {
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* ala  */
    MOL_RESIDUE_PROP_BASIC,          /* arg  */
    MOL_RESIDUE_PROP_POLAR,          /* asn  */
    MOL_RESIDUE_PROP_ACIDIC,         /* asp  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* cys  */
    MOL_RESIDUE_PROP_ACIDIC,         /* glu  */
    MOL_RESIDUE_PROP_POLAR,          /* gln  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* gly  */
    MOL_RESIDUE_PROP_BASIC,          /* his  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* ile  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* leu  */
    MOL_RESIDUE_PROP_BASIC,          /* lys  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* met  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* phe  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* pro  */
    MOL_RESIDUE_PROP_POLAR,          /* ser  */
    MOL_RESIDUE_PROP_POLAR,          /* thr  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* trp  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* tyr  */
    MOL_RESIDUE_PROP_HYDROPHOBIC,    /* val  */
    MOL_RESIDUE_PROP_UNKNOWN,        /* sol  */
    MOL_RESIDUE_PROP_UNKNOWN,        /* aden */
    MOL_RESIDUE_PROP_UNKNOWN,        /* guan */
    MOL_RESIDUE_PROP_UNKNOWN,        /* cyto */
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN,
    MOL_RESIDUE_PROP_UNKNOWN};

typedef enum MolAtomElementType {
  MOL_ATOM_ELEMENT_UNKNOWN,
  MOL_ATOM_ELEMENT_CARBON,
  MOL_ATOM_ELEMENT_HYDROGEN,
  MOL_ATOM_ELEMENT_NITROGEN,
  MOL_ATOM_ELEMENT_OXYGEN,
  MOL_ATOM_ELEMENT_SULFUR,
  MOL_ATOM_ELEMENT_PHOSPHORUS,
  MOL_ATOM_ELEMENT_IRON,
  MOL_ATOM_ELEMENT_ALL,
  MOL_ATOM_ELEMENT_MAX
  } MolAtomElementType;


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  MolAtom  ****                  *
 *                                                            *
 *------------------------------------------------------------*/

typedef struct MolAtom {
  int id;
  int orig_id;
  int res_seq;
  char insertion_code;
  int num;
  char *name;
  MolAtomElementType type;
  MolPoint3 pos;
  float temp;
  int het;
  char chain_id;
  MolResidueType res_type;
  MolResiduePropType res_prop;
  } MolAtom;

typedef struct MolAtomList {
  int num;
  int *list;
  } MolAtomList;

typedef struct MolAtomicCoords {
  int num;
  MolPoint3 *values;
  } MolAtomicCoords;

typedef struct MolAtomicStrucure {
  int num_atoms;
  int size;
  MolAtom *atoms;
  int atoms_grouped;
  MolAtomList *atom_groups[MOL_ATOM_ELEMENT_MAX];
  MolAtomicCoords coords;
  } MolAtomicStrucure;


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  MolResidue  ****               *
 *                                                            *
 *------------------------------------------------------------*/

typedef char (MolResName)[4];

typedef struct MolResidue {
  int id;
  char insertion_code;
  MolResidueType type;
  MolResiduePropType prop;
  MolAtomList atoms;
  int active;
  } MolResidue;

typedef struct MolResidueSeq {
  int num;
  int begin, end, index;
  MolResName *names;
  MolResidue **list;
  } MolResidueSeq;



/*------------------------------------------------------------*
 *                                                            *
 *                   **** PmHelixParams  ****                 *
 *                                                            *
 * monomer_dist: aka rise. dist btw monomers.                 *
 *------------------------------------------------------------*/

typedef struct MolHelixParams {
  MolPoint3 init_pt;
  MolPoint3 proj_pt;
  MolPoint3 proj_line[2];
  MolPoint3 center;
  float radius;
  MolVector3 axis, up, right;
  float ang0;
  int num_res;
  float monomer_dist;
  float res_per_turn;
  char *seq;
  int num_patterns;
  char patterns[30][80];
  int axis_rev;
  } MolHelixParams;

typedef struct MolHelix {
  int num;
  char id[3];
  char init_res_name[3];
  char init_chain_id;
  int init_seq_num;
  char term_res_name[3];
  char term_chain_id;
  int term_seq_num;
  int helixClass;
  int length;
  } MolHelix;

typedef struct MolSheet {
  int num;
  char id[3];
  int num_strands;
  char init_res_name[3];
  char init_chain_id;
  int init_seq_num;
  char init_icode;      
  char term_res_name[3];
  char term_chain_id;
  int term_seq_num;
  char end_icode;      
  int sense; 
  int length;
  } MolSheet;

typedef struct MolMoleculeProps {
  float center[3];
  float com[3];
  float inertia[3][3];
  float volume;
  float mass;
  float density;
  float area;
  float pscale[3], paxes[3][3];
  float pca_scale[3];
  MolPoint3 pca1, pca2, pca3;
  } MolMoleculeProps;


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  MolStructure  ****             *
 *                                                            *
 *------------------------------------------------------------*/

typedef struct MolChainTerm {
  int id;
  char chain_id;
  MolResName res_name;
  int res_seq;
  char insertion_code;
  struct MolChainTerm *next;
  } MolChainTerm;

typedef struct MolStructureSolvent {
  MolAtomicStrucure atomic_struct;
  int num_res;
  MolResidue *residues;
  int num_chains;
  char **chain_names;
  struct MolChain *chains;
  } MolStructureSolvent;

typedef struct MolStructure2ndHelix {
  int num;
  int size;
  MolHelix *list;
  } MolStructureSecHelix;

typedef struct MolStructure2ndSheet {
  int num;
  int size;
  MolSheet *list;
  } MolStructureSecSheet;

typedef struct MolStructureSec {
  MolStructureSecHelix helix;
  MolStructureSecSheet sheet;
  } MolStructureSec;

typedef struct MolStructure {
  char *name;
  int serial_number;
  MolMoleculeType mol_type;
  int num_chains;
  char **chain_names;
  struct MolChain *chains;
  MolChainTerm *terms;
  MolAtomicStrucure atomic_struct;
  MolAtomicStrucure het_atomic_struct;
  int extent_set;
  MolExtent extent;
  struct MolModel *model;
  int num_res;
  MolResidue *residues;
  MolStructureSec secondary;
  MolStructureSolvent solvent;
  struct MolStructure *next;
  } MolStructure;


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  MolChain  ****                 *
 *                                                            *
 *------------------------------------------------------------*/

typedef struct MolChain {
  char *name;
  char id;
  MolResidueSeq residues;
  int extent_set;
  MolExtent extent;
  MolStructure *struc;
  struct MolModel *model;
  MolMoleculeProps props;
  float color[3];
  int type;
  struct MolChain *next;
  } MolChain;


/*------------------------------------------------------------*
 *                                                            *
 *                       ****  MolModel  ****                 *
 *                                                            *
 *------------------------------------------------------------*/


typedef struct MolModel {
  const char *name;
  char *header;
  MolMoleculeType mol_type;
  MolExtent extent;
  int extent_set;
  int num_structures;
  MolStructure *structures;
  MolStructure *curr_structure;
  std::set<void*>* memory;
  } MolModel;



/*------------------------------------------------------------*
 *                       pdb parsing functions                *
 *------------------------------------------------------------*/

typedef struct MolDbProc {
  const char *name;
  void (*func)(FILE*, char*, MolModel*);
  } MolDbProc;

void 
mol_DbPdbAtomProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbCompProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbConectProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbHeaderProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbHelixProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbHetAtomProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbSeqresProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbSheetProc (FILE *fp, char *line, MolModel *model),
mol_DbPdbTerProc (FILE *fp, char *line, MolModel *model);

static MolDbProc  pdb_table[] = {
    {"atom",     mol_DbPdbAtomProc},
    {"compnd",   mol_DbPdbCompProc},
    {"conect",   mol_DbPdbConectProc},
    {"header",   mol_DbPdbHeaderProc},
    {"helix",    mol_DbPdbHelixProc},
    {"hetatm",   mol_DbPdbHetAtomProc},
    {"seqres",   mol_DbPdbSeqresProc},
    {"sheet",    mol_DbPdbSheetProc},
    {"ter",      mol_DbPdbTerProc},
    {NULL,       NULL}};


/*============================================================*
 *                  m e m o r y   a l l o c a t i o n         *
 *============================================================*/


//#define mem_Alloc(ptr, num, model) (ptr) = static_cast<typeof(ptr)>(mol_MemCallocFunc(#ptr, __FILE__, __LINE__,\
//                                                   (num), sizeof(*(ptr)), model))
//
//#define mem_Realloc(ptr,num,model) (ptr) = static_cast<typeof(ptr)>(mol_MemReallocFunc(#ptr, __FILE__, __LINE__, \
//                                                      (ptr), (num)*sizeof(*(ptr)), model))

// On windows, future standard keyword "typeof" is not recognized
#define mem_Alloc(ptr, num, model, type) (ptr) = static_cast<type>(mol_MemCallocFunc(#ptr, __FILE__, __LINE__,\
                                                   (num), sizeof(*(ptr)), model))

#define mem_Realloc(ptr, num, model, type) (ptr) = static_cast<type>(mol_MemReallocFunc(#ptr, __FILE__, __LINE__, \
                                                      (ptr), (num)*sizeof(*(ptr)), model))


/*============================================================*
 *                  p r o t o t y p e s                       *
 *============================================================*/

void
mol_AtomTypeConv (char *s, MolAtomElementType *type);

void
mol_ChainAtomsGet (MolChain *chain, int *p_num, int **p_list);

void
mol_ChainCreate (char chain_id, int num_res, int begin, int end, int rindex,
                 MolStructure *struc, int type, MolChain **p_chain);

void
mol_ChainResiduesGet (MolChain *chain, int *p_num, MolResidue ***p_list);

void 
mol_DbRead (const char *name, const char *file_name, MolDbType db_type, MolModel **p_model) ;

void 
mol_DbGromacsRead (MolModel *model, FILE *fp);

int
mol_DbLineGet (FILE *fp, char *line);

void
mol_DbLineStrip (char *s);

void
mol_DbRecordIntGet (char *line, int start, int end, int *val);

void
mol_DbRecordRealGet (char *line, int start, int end, float *val);

void
mol_DbRecordStrGet (char *line, int start, int end, char *str);

void
mol_DbMolSurfParse (FILE *fp, MolSurface *surf);

void
mol_DbMolSurfBinRead (FILE *fp, MolSurface *surf);

void *
mol_MemCallocFunc (const char *name, const char *file, int line, unsigned nelem, unsigned elsize, MolModel* model);

void *
mol_MemReallocFunc (const char *name, const char *file, int line, void *ptr, unsigned size, MolModel* model);

void
mol_MemFreeModel(MolModel* model);

void
mol_MolModelCreate (const char *name, MolModel **p_model);

void
mol_MolModelCurrStrucGet (MolModel *mobj, MolStructure **struc);

void
mol_MolModelHeaderSet (MolModel *model, char *hdr);

void
mol_MolModelStructureAdd (MolModel *model, MolStructure *struc);

void
mol_MolModelStructuresGet (MolModel *model, int *num, MolStructure **strucs);

void
mol_msg (const char *format, ...);

void
mol_ResTypeConv (char *s, MolResidueType *type);

void
mol_StructureAtomAdd (MolStructure *struc, int het, MolAtom *atom);

void
mol_StructureAtomsGet (MolStructure *struc, int *p_num, MolAtom **p_atoms);

void
mol_StructureChainsGet (MolStructure *struc, int *num, MolChain **chains);

void
mol_StructureCreate (char *name, MolStructure **p_struc, MolModel *model);

void
mol_StructureHelicesGet (MolStructure *struc, int *num, MolHelix **list);

void
mol_StructureModelGet (MolStructure *struc, MolModel **model);

void
mol_StructureResiduesGet (MolStructure *struc, int type, int *num_res, MolResidue **res);

void
mol_StructureSheetAdd (MolStructure *struc, MolSheet *sheet);

void
mol_StructureSheetsGet (MolStructure *struc, int *num, MolSheet **list);

void
mol_StructureTermAdd (MolStructure *struc, MolChainTerm *term, MolModel *model);

