
#include <iostream>

int readMRC(const char *fsp,int nodata, int n) { FILE *in;
//int EMData::readMRC(const char *fsp,int nodata, int n) { FILE *in;
	int i,j,l,pipe=0; int ord=0;	// 1 if the file order is MSB first int mord=0; // 1 if the machine order is MSB first float f,a,p; char s[800],*m2,w; unsigned char *cdata,*cdata2,c,t; short *sdata = 0;; unsigned short *usdata = 0; static char *str[] = { "char", "short", "float", "short complex", "float complex" }; static char *str2[] = { "" , "byte reversed" , "" }; static char *str3[] = { "" , "compressed" };
	flags=EMDATA_CHANGED; if (!mrch) mrch=(mrcH *)malloc(sizeof(mrcH));
	if (fsp==NULL) { in=stdin; } else {
	if (strcmp(&fsp[strlen(fsp)-3],".gz")==0 || strcmp(&fsp[strlen(fsp)-2],".Z")==0) { sprintf(s,"zcat %s",fsp);
	in=popen(s,"rb"); pipe=1;
	} if (strcmp(&fsp[strlen(fsp)-3],".bz")==0 || strcmp(&fsp[strlen(fsp)-4],".bz2")==0) {
	sprintf(s,"bzcat %s",fsp); in=popen(s,"rb"); pipe=1;
	} else in=fopen(fsp,"rb"); if (in==NULL) {
	error("Cannot open file",ERR_ERROR); zero(); return -1;
	} setPath(fsp);
	} if (!fread(mrch,sizeof(mrcH),1,in)) { setSize(10,10,1); zero(); return -1; }
	mrch->labels[0][79]=0; setName(mrch->labels[0]); if (name[0]=='!' && name[1]=='-') {
	i=sscanf(name+2," %f %f %f %f %f %f %f %f %f %f %f", ctf,ctf+1,ctf+2,ctf+3,ctf+4,ctf+5,ctf+6,ctf+7,ctf+8,ctf+9,ctf+10);
	if (i==11) flags|=EMDATA_HASCTF; else printf("Incomplete CTF info: %s\n",name);
	} if (name[0]=='!' && name[1]=='$') {
	sscanf(name+2," %f,%f",ctf,ctf+10);
	}
	// This swaps the byte order in the header on appropriate machines if (ByteOrder::is_machine_big_endian()) {
	mord=1;
	} else {
	mord=0;
	}
	m2=(char *)mrch; if (m2[0]==0 && m2[1]==0) ord=1; if (mord ^ ord) {
	}
	for (i=0; i<56*4; i+=4) { w=m2[i]; m2[i]=m2[i+3]; m2[i+3]=w; w=m2[i+1]; m2[i+1]=m2[i+2]; m2[i+2]=w;
	}
	if (mrch->ny>0xffff || mrch->nz>0xffff || mrch->mode>4 || mrch->nlabl>11) { setSize(10,10,1); printf("Invalid MRC file\n"); zero(); return -1;
	}
	int startz; if (mrch->mode==MODE_float_COMPLEX || mrch->mode==MODE_short_COMPLEX) mrch->nx*=2; nx=mrch->nx; ny=mrch->ny; if (n>=0 && n< mrch->nz) {
	} else {
	}
	nz=1; startz=n;
	if(n!=-1) printf("WARNING from readMRC: section n=%d is out of valid range [0,%d] and thus ignored, whole map will be read in\n",n,mrch->nz-1);
	nz=mrch->nz; startz=0;
	if (mrch->nz==1) sprintf(s,"Read a %s %s %dx%d 2D %s MRC File\n", str2[ord+mord],str3[pipe],nx,ny,str[mrch->mode]);
	else sprintf(s,"Read a %dx%dx%d 3D MRC File\n",nx,ny,nz);
	if (mrch->xlen==0) mrch->xlen=1.0; if (mrch->ylen==0) mrch->ylen=1.0; if (mrch->zlen==0) mrch->zlen=1.0;
	if (nodata) flags|=EMDATA_NODATA;
	// This reads 2d&3d arrays setSize(nx,ny,nz);
	if (!nodata) { portable_fseek(in,sizeof(mrcH)+mrch->nsymbt,SEEK_SET);
	getData(); cdata=(unsigned char *)rdata; sdata=(short *)rdata; usdata=(unsigned short *)rdata;
	switch(mrch->mode) { case MODE_char:
	portable_fseek(in,((EMAN_off_t)startz)*nx*ny,SEEK_CUR); flags&=~EMDATA_COMPLEX; for (j=0; j<ny; j++) {
	if (updProgress && updProgress(j,ny,"Reading MRC file (8 bit)")) break; if ((i=fread(&cdata[j*nx*nz],nx,nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_WARN);
	}
	} for (i=nx*ny*nz-1; i>=0; i--) rdata[i]=(float)cdata[i]/100.0-1.28; if (updProgress) updProgress(0,0,NULL); break;
	case MODE_short: portable_fseek(in,((EMAN_off_t)startz)*nx*ny*sizeof(short),SEEK_CUR); flags&=~EMDATA_COMPLEX;
	for (j=0; j<ny; j++) { if (updProgress && updProgress(j,ny,"Reading MRC file (16 bit)")) break; if ((i=fread(&sdata[j*nx*nz*sizeof(short)],nx*sizeof(short),nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_WARN); //printf ("ferror is %d", ferror(in));
	}
	}
	update(); if(fsp!=NULL){
	}
	} if ((ord+mord)&1) {
	} else {
	}
	for (i=nx*ny*nz-1; i>=0; i--) { t=cdata[i*2];
	cdata[i*2]=cdata[i*2+1]; cdata[i*2+1]=t; rdata[i]=(float)(sdata[i]);
	}
	for (i=nx*ny*nz-1; i>=0; i--) { rdata[i]=(float)(sdata[i]);
	}
	}
	} if (ord) {
	} else {
	for (i=nx*ny*nz-1; i>=0; i--) { rdata[i]=(float)((cdata[i*2]<<8)+cdata[i*2+1]);
	}
	for (i=nx*ny*nz-1; i>=0; i--) { rdata[i]=(float)((cdata[i*2+1]<<8)+cdata[i*2]);
	}
	//printf ("feof is %d", feof(in));
	if (updProgress) updProgress(0,0,NULL);
	break;
	case MODE_ushort: portable_fseek(in,((EMAN_off_t)startz)*nx*ny*sizeof(unsigned short),SEEK_CUR); flags&=~EMDATA_COMPLEX;
	for (j=0; j<ny; j++) { if (updProgress && updProgress(j,ny,"Reading MRC file (16 bit)")) break; if ((i=fread(&usdata[j*nx*nz*sizeof(unsigned short)],nx*sizeof(short),nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_WARN); //printf ("ferror is %d", ferror(in)); //printf ("feof is %d", feof(in));
	}
	} if ((ord+mord)&1) {
	} else {
	}
	for (i=nx*ny*nz-1; i>=0; i--) { t=cdata[i*2];
	cdata[i*2]=cdata[i*2+1]; cdata[i*2+1]=t; rdata[i]=(float)(usdata[i]);
	}
	for (i=nx*ny*nz-1; i>=0; i--) { rdata[i]=(float)(usdata[i]);
	}
	if (updProgress) updProgress(0,0,NULL);
	break;
	case MODE_float: {
	EMAN_off_t cur_offset = ((EMAN_off_t)startz)*nx*ny*sizeof(float); portable_fseek(in,cur_offset,SEEK_CUR); cdata=(unsigned char *)rdata; flags&=~EMDATA_COMPLEX;
	for (j=0; j<ny; j++) { if (updProgress && updProgress(j,ny,"Reading MRC file (float)")) break; if ((i=fread(&rdata[j*nx*nz],nx*sizeof(float),nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_WARN);
	}
	} // swap order if machine and file have different order if (ord+mord==1) {
	for (i=nx*ny*nz*4-4; i>=0; i-=4) { c=cdata[i]; cdata[i]=cdata[i+3]; cdata[i+3]=c; c=cdata[i+1]; cdata[i+1]=cdata[i+2]; cdata[i+2]=c;
	}
	} if (updProgress) updProgress(0,0,NULL); break;
	} case MODE_short_COMPLEX:
	portable_fseek(in,((EMAN_off_t)startz)*nx*ny*sizeof(short),SEEK_CUR); cdata=(unsigned char *)rdata; flags|=EMDATA_COMPLEX|EMDATA_RI; for (j=0; j<ny; j++) {
	if (updProgress && updProgress(j,ny,"Reading MRC file (16 bit complex)")) break; if ((i=fread(&cdata[j*nx*nz*2],nx*sizeof(short),nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_WARN);
	} if (updProgress) updProgress(0,0,NULL); break;
	case MODE_float_COMPLEX: portable_fseek(in,((EMAN_off_t)startz)*nx*ny*sizeof(float),SEEK_CUR); cdata=(unsigned char *)rdata; flags|=EMDATA_COMPLEX|EMDATA_RI; clearerr(in); for (j=0; j<ny; j++) {
	if (updProgress && updProgress(j,ny,"Reading MRC file (float complex)")) break; if ((i=fread(&rdata[j*nx*nz],nx*sizeof(float),nz,in))!=nz) {
	sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz); error(s,ERR_ERROR);
	}
	} if (ord+mord==1) {
	for (i=nx*ny*nz*4-4; i>=0; i-=4) { c=cdata[i]; cdata[i]=cdata[i+3]; cdata[i+3]=c; c=cdata[i+1]; cdata[i+1]=cdata[i+2]; cdata[i+2]=c;
	}
	} if (updProgress) updProgress(0,0,NULL); if (ferror(in)) { zero(); return -1; } break;
	default: printf("Bad MRC file\n");
	zero(); return -1;
	if (pipe) pclose(in); else fclose(in);
	}
	//setup pixel if possible pixel=mrch->xlen/(float)(nx); xorigin=mrch->xorigin; yorigin=mrch->yorigin; zorigin=mrch->zorigin; parent=NULL;
	vol_table = 0; pathnum=0;
	// somebody changed the writeMRC to use nx instead of nx-1, so this is update to use nx too for self consistence
	dx=dy=dz=0; setRAlign(0,0,0); atime=time(0); doneData(); update();
	// this is to filp the phase for conforming to MRC convention if (isComplex()) {
	ap2ri(); getData(); for(i=0;i<nx*ny*nz;i+=2)rdata[i+1]*=-1; doneData();
	} return 0;
}
