package var;

import dna.Gene;


public class VarLine extends Variation{
	
	public static final long serialVersionUID = -4089933371294357462L;
	
//	>locus  ploidy  haplotype       chromosome      begin   end     varType reference       alleleSeq       totalScore      hapLink xRef

	public VarLine(){}

	public VarLine(String s, float version){
		String[] line=s.split("\t", -1);
		
		for(int i=0; i<line.length; i++){
			line[i]=line[i].trim();
			if(line[i].length()<1){
				line[i]=null;
			}
		}
		
		
		
//		varType=(byte)find((line.length>6 ? line[6] : "null"), varTypeMap);
		Byte b=varTypeMap2.get(line.length>6 ? line[6] : "null");
		assert(b!=null) : "Can't find "+line[6]+" in "+varTypeMap2.keySet()+"\n\nLine: "+s+"\n";
		varType=b;
		
		
//		locus=Integer.parseInt(line[0]);
		
		b=(Byte)ploidyMap.get(line[1]);
		assert(b!=null) : "\n\n"+line[1]+"\n\n"+s+"\n\n";
		ploidy=b;
		
		
		haplotype=(byte)find(line[2], haploMap);
		assert(haplotype>=0) : line[2];
		
		chromosome=Gene.toChromosome(line[3]);
		assert(chromosome>0) : line[3]+" -> "+line[3].substring(3);
		
		beginLoc=Integer.parseInt(line[4]);
		int tempInt=Integer.parseInt(line[5])-1; //Note:  0,1 based
		tempInt=max(tempInt, beginLoc);
		endLoc=tempInt;
		
		String temp;
		
		temp=line.length>7 ? line[7] : null;
		if("?".equals(temp)){temp=null;}
		ref=temp;
		
		temp=line.length>8 ? line[8] : null;
		if("?".equals(temp)){temp=null;}
		call=temp;
		
		
		if(version<2){
			
			totalScore=((line.length<=9 || line[9]==null || line[9].length()<1) ? -1 : Integer.parseInt(line[9]));
			hapLink=((line.length<=10 || line[10]==null || line[10].length()<1) ? -1 : Integer.parseInt(line[10]));

			assert(beginLoc<=endLoc) : s;
			
			//		System.out.println("\n"+this+"\n"+new Variation(this)+"\n");
		}else{

//			return "#locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\t
//			varScoreVAF\tvarScoreEAF\tvarQuality\thapLink\txRef
			
			int varScoreVAF=((line.length<=9 || line[9]==null || line[9].length()<1) ? -1 : Integer.parseInt(line[9]));
			int varScoreEAF=((line.length<=10 || line[10]==null || line[10].length()<1) ? -1 : Integer.parseInt(line[10]));
			byte VQ=((line.length<=11 || line[11]==null || line[11].length()<1) ? (byte)0 : (byte)find(line[11], VQARRAY));
			
			totalScore=varScoreVAF;
			hapLink=((line.length<=12 || line[12]==null || line[12].length()<1) ? -1 : Integer.parseInt(line[12]));
			
			assert(beginLoc<=endLoc) : s;
			
//			System.out.println("\n"+this+"\n"+new Variation(this)+"\n");
		}
		
		assert(!((varType==Variation.INS || varType==Variation.DELINS || varType==Variation.SNP)
				&& call==null)) : "\nversion="+version+"\n"+s+"\n"+line+"\nline.ref="+ref+"\nline.call="+call+"\nref="+ref+"\ncall="+call;
		
		intern();
	}
	
	@Override
	public VarLine clone(){
		VarLine v=null;
//		try {
//			v=(VarLine) super.clone();
//		} catch (CloneNotSupportedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		v=(VarLine) super.clone();
		return v;
	}
	
	public VarLine[] splitLine(){
		assert(haplotype==3) : this;
		VarLine[] r=new VarLine[2];
		r[0]=this.clone();
		r[1]=this.clone();
		assert(this.equals(r[0]) && r[0].equals(this));
		r[0].haplotype=1;
		r[1].haplotype=2;
		return r;
	}
	
	public VarLine spawnEqualPoint(){
		assert(this.isPoint());
		VarLine v=this.clone();
		v.varType=REFPOINT;
		v.call=v.ref=null;
		return v;
	}
	
	public static VarLine makeEqualPoint(byte chrom, int loc, byte hap){
		VarLine v=new VarLine();
		v.chromosome=chrom;
		v.beginLoc=loc;
		v.endLoc=loc;
		v.haplotype=hap;
		v.varType=REFPOINT;
		return v;
	}
	
	@Override
	public String toSuperString(){return super.toString();}
	
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder(256);

//		sb.append(locus+"\t");
		sb.append(ploidyMap.get(ploidy)+"\t");
		sb.append(haploMap[haplotype]+"\t");
		sb.append("chr"+Gene.chromCodes[chromosome]+"\t");
		sb.append(beginLoc+"\t");
		sb.append(endLoc+"\t");
		
		sb.append(varTypeMap[varType]+"\t");
		sb.append((ref==null ? "" : ref)+"\t");
		sb.append((call==null ? "" : call)+"\t");
		sb.append((totalScore==-1 ? "" : totalScore)+"\t"); //TODO: Note the collision with a true -1
		sb.append((hapLink==-1 ? "" : hapLink)+"\t"); //TODO "
		
		return sb.toString();
	}
	
	public static String sourceHeader(){
		return "#locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\ttotalScore\thapLink\txRef";
//		return "#locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\t
//			varScoreVAF\tvarScoreEAF\tvarQuality\thapLink\txRef
	}
	//locus	ploidy	allele	chromosome	begin	end	varType	reference	alleleSeq	varScoreVAF	varScoreEAF	varQuality	hapLink	xRef
	
	@Override
	public String toSourceString(){
		StringBuilder sb=new StringBuilder(256);

		sb.append(0+"\t");
		sb.append(ploidyMap.get(ploidy)+"\t");
		sb.append(haploMap[haplotype]+"\t");
		sb.append("chr"+Gene.chromCodes[chromosome]+"\t");
		sb.append(beginLoc+"\t");
		
		if(varType==INS){
			sb.append(beginLoc+"\t");
		}else{
			sb.append((endLoc+1)+"\t");
		}
		
		sb.append(varTypeMap[varType]+"\t");
		sb.append((ref==null ? "" : ref)+"\t");
		sb.append((call==null ? "" : call)+"\t");
		sb.append((totalScore==-1 ? "" : totalScore)+"\t"); //TODO: Note the collision with a true -1
		sb.append((hapLink==-1 ? "" : hapLink)+"\t"); //TODO "
		
		return sb.toString();
	}
	
	
	@Override
	public String toShortString(){
		StringBuilder sb=new StringBuilder(256);
		
		sb.append(haploMap[haplotype]);
		while(sb.length()<3){sb.append(' ');}
		sb.append('\t');
		sb.append(locationString()+"\t");
		
		sb.append(varTypeMap[varType]+"\t");
		sb.append((ref==null ? "" : ref)+"\t");
		sb.append((call==null ? "" : call)+"\t");
//		sb.append((totalScore==-1 ? "" : totalScore)+"\t"); //TODO: Note the collision with a true -1
//		sb.append((hapLink==-1 ? "" : hapLink+"\t")); //TODO "
		
		return sb.toString();
	}
	
	@SuppressWarnings("unused")
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	
	@Override
	public int compareTo(Variation other) {
		if(other.getClass()==VarLine.class){
			return compareTo((VarLine)other);
		}
		return super.compareTo(other);
	}
	
	public int compareTo(VarLine other) {
		int x=super.compareTo((Variation)other);
		if(x!=0){return x;}
		return haplotype-other.haplotype;
	}
	
	@Override
	public boolean equals(Object other){
		if(other.getClass()==VarLine.class){
			return equals((VarLine)other);
		}
		return super.equals(other);
	}
	
	public boolean equals(VarLine other){
		return compareTo(other)==0;
	}
	
	@Override
	public boolean equals(Variation other){
		return super.equals(other);
	}
	
	public byte ploidy;

	/** Which copy this is on */
	public byte haplotype;
	public int totalScore;
	public int hapLink;
	
	public static final String[] VQARRAY=new String[] {"?", "VQLOW", "VQHIGH"};
	
}
