package dna;
import java.io.Serializable;

import shared.Shared;
import shared.Tools;


public class Gene implements Comparable<Gene>, Serializable{
	
//	/**
//	 *
//	 */
	private static final long serialVersionUID = -1342555621377050981L;
	
	
	public Gene(){
		chromosome=-1;
//		nc_accession=null;
		symbol=null;
		proteinAcc=null;
		id=-1;
		mrnaAcc=null;
		status=-1;
		completeness=-1;
		strand=-1;
		codeStart=txStart=-1;
		codeStop=txStop=-1;
		exons=null;
		cdsStartStat=-1;
		cdsEndStat=-1;
		exonFrames=null;
		txLength=-1;
		codeLength=-1;
		exonLength=-1;
		exonCodeLength=-1;
		aaLength=-1;
		utrLength5prime=-1;
		utrLength3prime=-1;
		readCorrectly=false;
		untranslated=false;
		pseudo=false;
		description=null;
		fullDescription=null;
		valid=true;
		primarySource=-1;
	}
	
	public Gene(byte chrom, byte strand_, int txStart_, int txStop_, int cdStart_, int cdStop_, int gid,
			String name_, String trans_, String protTrans_, String status_, String completeness_,
			Exon[] exons_, boolean untran, boolean pseudo_, boolean valid_,
			String primarySource_, String descript_, String fullDescript_){
		
		chromosome=chrom;
//		nc_accession=null;
		symbol=name_;
		id=gid;
		mrnaAcc=((trans_==null || trans_.length()<1 || trans_.equals("-")) ? null : trans_);
		proteinAcc=((protTrans_==null || protTrans_.length()<1 || protTrans_.equals("-")) ? null : protTrans_);
		
		primarySource=primarySource_==null ? -1 : (byte)Tools.find(primarySource_, sourceCodes);
		description=descript_;
		fullDescription=fullDescript_;
		
		
		status=status_==null ? -1 : (byte)Tools.find(status_, statusCodes);
		completeness=completeness_==null ? -1 : (byte)Tools.find(completeness_, completenessCodes);
		strand=strand_;
		
		exons=exons_;

		txStart=txStart_;
		txStop=txStop_; //Assuming pure 0-based numbering.
		codeStart=cdStart_;
		codeStop=cdStop_; //Assuming pure 0-based numbering.

		assert(codeStart>=txStart) : "("+txStart+", "+txStop+"), ("+codeStart+", "+codeStop+") for "+mrnaAcc;
		assert(codeStop<=txStop) : "("+txStart+", "+txStop+"), ("+codeStart+", "+codeStop+") for "+mrnaAcc;
		

//		cdsStartStat=(byte)find("?", endStatCodes);
//		cdsEndStat=(byte)find("?", endStatCodes);
		cdsStartStat=-1;
		cdsEndStat=-1;
		
		exonFrames=null;
		
		txLength=txStop-txStart+1;
		codeLength=(codeStop==codeStart ? 0 : codeStop-codeStart+1);
		
		untranslated=untran;
		pseudo=pseudo_;
		
		int eLen=0, ecLen=0, utr0=0, utr2=0;
		
		if(exons!=null){

			for(Exon e : exons){

				utr0+=max(0, min(e.b, codeStart)-e.a);
				utr2+=max(0, e.b-max(e.a, codeStop));
				
				int len=e.b-e.a+1;
				eLen+=len;
				len=(min(e.b, codeStop)-max(e.a, codeStart));
				len=max(0, len+1);
				ecLen+=len;
			}
		}
		
		
		exonLength=(eLen<2 ? 0 : eLen);
		exonCodeLength=(codeLength<1 || exonLength<1 ? 0 : ecLen);
		aaLength=exonCodeLength/3-1;

		assert(exonLength>=exonCodeLength) : exonLength+", "+codeLength+", "+exonCodeLength+"\n"+this+"\n";
		assert(codeLength>=exonCodeLength) : exonLength+", "+codeLength+", "+exonCodeLength+"\n"+this+"\n";
		
		//assert(exonCodeLength%3 == 0); //This should be true with a correct database
		
		if(strand==Shared.PLUS){
			utrLength5prime=untranslated ? 0 : utr0;
			utrLength3prime=untranslated ? 0 : utr2;
		}else{
			utrLength5prime=untranslated ? 0 : utr2;
			utrLength3prime=untranslated ? 0 : utr0;
		}

		//System.err.println(name+", "+exonLength+", "+exonCodeLength+(exons==null ? "" : ", "+exons.length));
		
		readCorrectly=true;
		valid=(readCorrectly && valid_);
	}
	
	
	public Gene merge(Gene g){

		assert((exons==null && g.exons==null) ||
				(exons!=null && g.exons!=null && exons.length==g.exons.length));
//		assert(exonLength==g.exonLength);
		assert(Math.abs(exonLength-g.exonLength)<=8) : "\n\n"+this+"\n\n"+g+"\n\n";
		assert(strand==g.strand);
//		assert(codeStart==g.codeStart);
//		assert(codeStop==g.codeStop);
		
		String Xsymbol=symbol;
		String XproteinAcc=proteinAcc;
		int Xid=id;
		String XmrnaAcc=mrnaAcc;
		int Xstatus=status;
		int Xcompleteness=completeness;
		int XcodeStart=codeStart;
		int XcodeStop=codeStop;
		int XtxStart=txStart;
		int XtxStop=txStop;
		int XcdsStartStat=cdsStartStat;
		int XcdsEndStat=cdsEndStat;
		byte[] XexonFrames=exonFrames;
		int XtxLength=txLength;
		int XcodeLength=codeLength;
		int XexonLength=exonLength;
		int XexonCodeLength=exonCodeLength;
		int XaaLength=aaLength;
		int XutrLength5prime=utrLength5prime;
		int XutrLength3prime=utrLength3prime;
//		boolean XreadCorrectly=readCorrectly;
		boolean Xuntranslated=untranslated;
		boolean Xpseudo=pseudo;
		String Xdescription=description;
		String XfullDescription=fullDescription;
		boolean Xvalid=valid;
		
		assert(untranslated || g.untranslated || g.codeStart>=txStart) : "\n"+this+"\n\n"+g;
		assert(untranslated || g.untranslated || g.codeStop<=txStop) : "\n"+this+"\n\n"+g;
		
		if(Xsymbol==null){Xsymbol=g.symbol;}
		if(XproteinAcc==null){XproteinAcc=g.proteinAcc;}
		if(Xid<0){Xid=g.id;}
		if(XmrnaAcc==null){XmrnaAcc=g.mrnaAcc;}
		if(Xstatus<0){Xstatus=g.status;}
		if(Xcompleteness<0){Xcompleteness=g.completeness;}
		
		
		if(XcodeStart==XcodeStop && g.codeStart<g.codeStop){
			assert(g.codeStart>=txStart);
			assert(g.codeStop<=txStop);
			XcodeStart=g.codeStart;
			XcodeStop=g.codeStop;
		}
		
		//These two should never happen...
		if(XtxStart<0){XtxStart=g.txStart;}
		if(XtxStop<0){XtxStop=g.txStop;}
		
		if(XcdsStartStat<0){XcdsStartStat=g.cdsStartStat;}
		if(XcdsEndStat<0){XcdsEndStat=g.cdsEndStat;}
		if(XexonFrames==null){XexonFrames=g.exonFrames;}
		if(XtxLength<0){XtxLength=g.txLength;}
		if(XcodeLength<0){XcodeLength=g.codeLength;}
		if(XexonLength<0){XexonLength=g.exonLength;}
		if(XexonCodeLength<0){XexonCodeLength=g.exonCodeLength;}
		if(XaaLength<0){XaaLength=g.aaLength;}
		if(XutrLength5prime<0){XutrLength5prime=g.utrLength5prime;}
		if(XutrLength3prime<0){XutrLength3prime=g.utrLength3prime;}
		if(Xdescription==null){Xdescription=g.description;}
		if(XfullDescription==null){XfullDescription=g.fullDescription;}
		
//		if(XreadCorrectly){}
//		if(Xuntranslated){}
//		if(Xpseudo){}
//		if(Xvalid){}
		
		//TODO Note that the readCorrectly field gets lost here
		Gene out=new Gene(chromosome, strand, XtxStart, XtxStop, XcodeStart, XcodeStop, Xid,
				symbol, XmrnaAcc, XproteinAcc,
				Xstatus< 0 ? null : statusCodes[Xstatus], Xcompleteness<0 ? null : completenessCodes[Xcompleteness],
				exons, Xuntranslated, Xpseudo, Xvalid, sourceCodes[primarySource], Xdescription, XfullDescription);
		
		return out;
	}
	
	public static byte toStrand(String s){
		assert(s!=null && s.length()==1);
		final char c=s.charAt(0);
		if(c=='+'){return Shared.PLUS;}
		else if(c=='-'){return Shared.MINUS;}
		else if(c=='?'){return 2;}
		throw new RuntimeException("Unknown strand: "+s);
	}
	
	public static int toChromosome(final String s){
////		assert(false) : s;
//		String s2=s;
//		if(s2.endsWith("random")){s2="U";}
//		if(s2.startsWith("chr")){s2=s2.substring(3);}
//		if(s2.equals("MT")){s2="M";}
////		int loc=find2(s2.toUpperCase(), chromCodes);
//		int loc=find3(s2.toUpperCase(), chromCodes);
//
//		if(loc<0){
//			if(!Tools.isDigit(s2.charAt(0))){
//				loc=find3("U", chromCodes);
//			}else{
//				try {
//					loc=Integer.parseInt(s2);
//				} catch (NumberFormatException e) {
//					throw new RuntimeException(e);
//				}
//				assert(loc>=23 && loc<=26) : loc+", "+s;
//			}
//		}
//		assert(loc>=0) : s;
//		return loc;
		
		String s2=s;
		if(s2.startsWith("chr")){s2=s2.substring(3);}
		int loc=Integer.parseInt(s2);
		
		assert(loc>=0) : s;
		return loc;
	}
	
	public static int toBuild(final String s){
		String s2=s.toLowerCase();
		if(s2.startsWith("build")){s2=s2.substring(5);}
		else if(s2.startsWith("b")){s2=s2.substring(1);}
		else if(s2.startsWith("hg")){s2=s2.substring(1);}
		
		if(s2.startsWith("=")){s2=s2.substring(1);}
		
		assert(Tools.isDigit(s2.charAt(0))) : s;
		
		return Integer.parseInt(s2);
	}
	
	private void fillExons(String eStarts, String eEnds, byte chr, byte str){
		String[] s1=eStarts.split(",");
		String[] s2=eEnds.split(",");
		
		int last=-1;
		
		for(int i=0; i<s1.length; i++){
			int a=Integer.parseInt(s1[i]);
			int b=Integer.parseInt(s2[i])-1; //Note the -1 for 0-based numbering.
			assert(a>last) : eStarts;
			last=a;
			
			boolean cds=overlap(a, b, codeStart, codeStop);
			boolean utr=(a<codeStart || b>codeStop);
			
			Exon key=new Exon(a, b, chr, str, utr, cds);
			Exon value=Exon.table.get(key);
			if(value==null){
				value=key;
				Exon.table.put(key, key);
			}
			exons[i]=value;
		}
	}
	
	private Exon[] fillExonsCCDS(String estring, byte chr, byte str){
		String[] intervals=estring.replace("[","").replace("]","").replace(" ","").split(",");
		
		int last=-1;
		
		Exon[] array=new Exon[intervals.length];
		
		for(int i=0; i<intervals.length; i++){
			String[] temp=intervals[i].split("-");
			int a=Integer.parseInt(temp[0]);
			int b=Integer.parseInt(temp[1]); //Note the pure 0-based numbering.
			assert(a>last) : estring;
			last=a;
			
			boolean cds=overlap(a, b, codeStart, codeStop);
			boolean utr=(a<codeStart || b>codeStop);
			
			Exon key=new Exon(a, b, chr, str, utr, cds);
			Exon value=Exon.table.get(key);
			if(value==null){
				value=key;
				Exon.table.put(key, key);
			}
			array[i]=value;
		}
		return array;
	}
	
	public int toGeneRelativeOffset(int index){
		
		int off=0;

		if(strand==Shared.PLUS){

			//		System.out.println();
			for(Exon e : exons){
				//			System.out.print(e+" * ");

				int temp=0;
				if(e.intersects(index)){
					temp=(int)(index-e.a);
				}else if(e.a>index){
					break;
				}else{
					temp=e.length();
				}
				assert(temp<=e.length()) : index +" \t "+e+" \t "+temp+" \t "+e.length();
				assert(temp>=0) : index+", "+e;
				off+=temp;
			}

		}else if(strand==Shared.MINUS){
			for(int i=exons.length-1; i>=0; i--){
				Exon e=exons[i];

				int temp=0;
				if(e.intersects(index)){
					temp=(int)(e.b-index);
				}else if(e.b<index){
					break;
				}else{
					temp=e.length();
				}
				assert(temp<=e.length()) : index +" \t "+e+" \t "+temp+" \t "+e.length();
				assert(temp>=0) : index+", "+e;
				off+=temp;
			}

		}else{assert false : strand;}
		
		return off;
	}
	
	public int[] toExonRelativeOffset(int index){
		
		int ex=0;
		int off=0;

		if(strand==0){

			//		System.out.println();
			for(Exon e : exons){
				//			System.out.print(e+" * ");

				int temp=0;
				if(e.intersects(index)){
					temp=(int)(index-e.a);
				}else if(e.a>index){
					break;
				}else{
					ex++;
				}
				assert(temp<=e.length()) : index +" \t "+e+" \t "+temp+" \t "+e.length();
				assert(temp>=0) : index+", "+e;
				off=temp;
			}

		}else if(strand==1){
			for(int i=exons.length-1; i>=0; i--){
				Exon e=exons[i];

				int temp=0;
				if(e.intersects(index)){
					temp=(int)(e.b-index);
				}else if(e.b<index){
					break;
				}else{
					ex++;
				}
				assert(temp<=e.length()) : index +" \t "+e+" \t "+temp+" \t "+e.length();
				assert(temp>=0) : index+", "+e;
				off=temp;
			}

		}else{assert false : strand;}
		
//		if((index-143053138)>-3 && (index-143053138)<3){
//			assert(false) : ("\n\nLooking for "+index+" in\n"+this+
//					"\n\nwith exons\n"+Arrays.toString(exons)+"\n\nResult: "+off+"\n\n");
//		}
//
//		if((index-143053111)>-10 && (index-143053111)<10){
//			assert(false) : ("\n\nLooking for "+index+" in\n"+this+
//					"\n\nwith exons\n"+Arrays.toString(exons)+"\n\nResult: "+off+"\n\n");
//		}
		
//		if(off==1 && exons[exons.length-1].b==143053111){
//			assert(false) : ("\n\nLooking for "+index+" in\n"+this+
//					"\n\nwith exons\n"+Arrays.toString(exons)+"\n\nResult: "+off+"\n\n");
//		}
		
		//		System.out.println();
		return new int[] {ex, off};
	}
	
	
	public boolean isHypothetical(){
		return isHypothetical(symbol);
	}
	
	
	public static boolean isHypothetical(String s){
		if(s==null){return false;}
		if(s.startsWith("C") && s.contains("orf")){return true;}
		if(s.length()>=4 && s.startsWith("LOC") && Tools.isDigit(s.charAt(3))){return true;}
		return false;
	}
	
	
	public boolean isNormalGene(){
		return valid && !untranslated && !pseudo && !isHypothetical();
	}
	

	public boolean intersectsTx(int point){
		return point>=txStart && point<=txStop;
	}
	public boolean intersectsTr(int point){
		assert(!untranslated);
		return (untranslated ? false : point>=translationStart() && point<=translationStop());
	}
	public boolean intersectsCode(int point){
//		assert(!untranslated) : "point = "+point+"\ngene = "+this;
//		return (untranslated ? false : point>=codeStart && point<=codeEnd);
		return (untranslated ? intersectsTx(point) : point>=codeStart && point<=codeStop);
	}
	public boolean intersectsExon(int point){
		for(Exon e : exons){
			if(e.intersects(point)){return true;}
		}
		return false;
	}
	
	/** Note that this skips code intersection checking for untranslated genes. */
	public boolean intersectsCodeAndExon(int point){
		if(!untranslated && !intersectsCode(point)){return false;}
		for(Exon e : exons){
			if(e.intersects(point)){return true;}
		}
		return false;
	}
	
	
	/** Note that this skips code intersection checking for untranslated genes. */
	public boolean intersectsCodeAndExon(int a, int b){
		if(!untranslated && !intersectsCode(a, b)){return false;}
		for(Exon e : exons){
			if(e.intersects(a, b)){return true;}
		}
		return false;
	}
	
	/** Note that this skips code intersection checking for untranslated genes. */
	public boolean intersectsIntron(int a, int b){
		if(exons==null || exons.length<2){return false;}
		if(!overlap(a, b, exons[0].a, exons[exons.length-1].b)){return false;}
		for(int i=1; i<exons.length; i++){
			Exon e1=exons[i-1];
			Exon e2=exons[i];
			assert(e1.b<e2.a) : "\n"+e1+"\n"+e2+"\n"+this+"\n";
			
			assert(a<=b && e1.b+1<=e2.a-1) : "\n"+e1+"\n"+e2+"\n"+this+"\n";
			
			if(overlap(a, b, e1.b+1, e2.a-1)){return true;}
		}
		return false;
	}
	
	/** Note that this skips code intersection checking for untranslated genes. */
	public boolean isDeepIntronic(int a, int b, int distFromEnds){
		if(exons==null){return false;}
		for(int i=1; i<exons.length; i++){
			Exon e1=exons[i-1];
			Exon e2=exons[i];
			assert(e1.b<e2.a) : "\n"+e1+"\n"+e2+"\n"+this+"\n";
			if(a>=e1.b+distFromEnds && b<=e2.a-distFromEnds){return true;}
		}
		return false;
	}
	
	public boolean intersectsSplice(int a, int b){
		assert(b>=a);
		if(exons==null || exons.length<2){return false;}
		if(b<txStart || a>txStop){return false;}
		for(Exon e : exons){
			if(e.a>=a && e.a<=b){return true;}
			if(e.b>=a && e.b<=b){return true;}
		}
		return false;
	}
	
	public boolean intersectsNearby(int a, int b){
		return intersectsCodeAndExon(a-NEAR, b+NEAR);
	}
	
	private static int closestToPoint(int a, int b, int point){
		int a2=(a>point ? a-point : point-a);
		int b2=(b>point ? b-point : point-b);
		return a2<b2 ? a : b;
	}
	
	/**
	 * @param a
	 * @param b
	 * @return {
	 * 	distance,<br>
	 * 	nearest exon number (-1 means coding start or stop),<br>
	 * 	side (0 means start, 1 means stop),<br>
	 * 	position (1 means inside, 2 means outside, 3 means both),<br>
	 *  site coordinate
	 *  }
	 */
	public int[] nearestSpliceSite(int a, int b){
		
		int bestDist=999999999;
		int nearestExon=-1;
		int side=-1;
		int position=0;
		int bestSite=-1;
		
		boolean strictlyIntronic=this.isDeepIntronic(a, b, 1);
		
		if(!strictlyIntronic){
			{
				int point=codeStart;
				int x=Exon.distToPoint(a, b, point);
				if(x<bestDist){
					bestDist=x;
					bestSite=point;
					nearestExon=-1;
					position=0;
					if(a<point){position|=2;}
					if(b>=point){position|=1;}
					side=(strand==Shared.PLUS ? 0 : 1);
					if(strand==Shared.PLUS){
						side=0;
					}else if(strand==Shared.MINUS){
						side=1;
					}
				}

				point=codeStop;
				x=Exon.distToPoint(a, b, point);
				if(x<bestDist){
					bestDist=x;
					bestSite=point;
					nearestExon=-1;
					position=0;
					if(b>point){position|=2;}
					if(a<=point){position|=1;}
					side=(strand==Shared.PLUS ? 1 : 0);
				}
			}
		}
		
		for(int i=0; i<exons.length; i++){
			Exon e=exons[i];
			
			int point=e.a;
			int x=Exon.distToPoint(a, b, point);
			if(x<bestDist){
				bestDist=x;
				bestSite=point;
				nearestExon=i;
				side=(strand==Shared.PLUS ? 0 : 1);
				position=0;
				if(a<point){position|=2;}
				if(b>=point){position|=1;}
			}
			
			point=e.b;
			x=Exon.distToPoint(a, b, point);
			if(x<bestDist){
				bestDist=x;
				bestSite=point;
				nearestExon=i;
				side=(strand==Shared.PLUS ? 1 : 0);
				position=0;
				if(b>point){position|=2;}
				if(a<=point){position|=1;}
			}
		}
		
		if(nearestExon>=0 && strand==Shared.MINUS){
			nearestExon=exons.length-nearestExon-1;
		}
		
		return new int[] {bestDist, nearestExon, side, position, bestSite};
	}
	
	
	
	public boolean intersectsTx(int a, int b){
		assert(a<=b);
		return overlap(a, b, txStart, txStop);
	}
	public boolean intersectsTr(int a, int b){
		assert(a<=b);
		assert(!untranslated);
		return (untranslated ? false : overlap(a, b, translationStart(), translationStop()));
	}
	public boolean intersectsCode(int a, int b){
		assert(a<=b);
//		assert(!untranslated) : "a="+a+", b="+b+"\ngene = "+this;
//		return (untranslated ? false : overlap(a, b, codeStart, codeEnd));
		return (untranslated ? intersectsTx(a, b) : overlap(a, b, codeStart, codeStop));
	}
	public boolean intersectsExon(int a, int b){
//		if(!intersectsCode(a, b)){return false;}
		assert(a<=b);
		for(Exon e : exons){
			if(e.intersects(a, b)){return true;}
		}
		return false;
	}
	public boolean intersectsUTR(int a, int b){
		if(!intersectsTx(a,b)){return false;}
		if(untranslated){return true;}
		if(overlap(a, b, txStart, codeStart)){return true;}
		if(overlap(a, b, codeStop, txStop)){return true;}
		return false;
	}
	/** Downstream */
	public boolean intersectsUTR3(int a, int b){
		if(!intersectsTx(a,b)){return false;}
		if(untranslated){return false;}
		if(strand==Shared.MINUS){
			if(overlap(a, b, txStart, codeStart)){return true;}
		}else{
			if(overlap(a, b, codeStop, txStop)){return true;}
		}
		return false;
	}
	/** Upstream */
	public boolean intersectsUTR5(int a, int b){
		if(!intersectsTx(a,b)){return false;}
		if(untranslated){return false;}
		if(strand==Shared.PLUS){
			if(overlap(a, b, txStart, codeStart)){return true;}
		}else{
			if(overlap(a, b, codeStop, txStop)){return true;}
		}
		return false;
	}
	
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	public static final String header(){
		return "#chrom\tsymbol\tgeneId\tmrnaAcc\tproteinAcc" +
				"\tstrand\tcodeStart\tcodeStop\ttxStart\ttxStop" +
				"\t(UNTRANSLATED?)\t(PSEUDOGENE?)\tstatus\tcompleteness\tsource" +
				"\t[exon0start-exon0stop, ...exonNstart-exonNstop]" +
				"\tfullName\tdescription";
	}
	
//	public CharSequence toRefSeqFormat(){
//		return driver.ToRefGeneFormat.format(this);
//	}
	
	@Override
	public String toString(){
		
		StringBuilder sb=new StringBuilder(256);

		sb.append(chromosome+"\t");
		sb.append(symbol+"\t");
		sb.append(id+"\t");
		sb.append(mrnaAcc+"\t");
		assert(proteinAcc==null || !proteinAcc.equals("null"));
		sb.append((proteinAcc==null ? "" : proteinAcc)+"\t");
		sb.append(Shared.strandCodes[strand]+"\t");
		sb.append(codeStart+"\t");
		sb.append(codeStop+"\t");
		sb.append(txStart+"\t");
		sb.append(txStop);
		
		sb.append("\t"+(untranslated ? "UNTRANS" : ""));
		sb.append("\t"+(pseudo ? "PSEUDO" : ""));
		
		sb.append("\t"+(status>=0 ? statusCodes[status] : ""));
		sb.append("\t"+(completeness>=0 ? completenessCodes[completeness] : ""));
		sb.append("\t"+(primarySource>=0 ? sourceCodes[primarySource] : ""));
		
		sb.append("\t[");
		String comma="";
		for(int i=0; exons!=null && i<exons.length; i++){
			sb.append(comma+exons[i].a+"-"+exons[i].b);
			comma=", ";
		}
		sb.append("]");

		assert(description==null || (!description.equals("null") && description.length()>0));
		sb.append("\t"+(description==null ? "" : description));
		
		assert(fullDescription==null || (!fullDescription.equals("null") && fullDescription.length()>0));
		sb.append("\t"+(fullDescription==null ? "" : fullDescription));
		
		String s=sb.toString();
		return Character.isWhitespace(s.charAt(0)) ? s : s.trim();
	}
	
	public String toShortString(){
		
		StringBuilder sb=new StringBuilder(256);

		sb.append("chr"+chromosome+"\t");
		sb.append(symbol+"\t");
		sb.append(mrnaAcc+"\t");
		sb.append(Shared.strandCodes[strand]+"\t");
		sb.append("("+codeStart+" - "+codeStop+")");
		return sb.toString();
	}
	
	@Override
	public int compareTo(Gene other){
		if(chromosome<other.chromosome){return -1;}
		if(chromosome>other.chromosome){return 1;}
		
		if(txStart<other.txStart){return -1;}
		if(txStart>other.txStart){return 1;}
		
		if(txStop<other.txStop){return -1;}
		if(txStop>other.txStop){return 1;}
		
		if(codeStart<other.codeStart){return -1;}
		if(codeStart>other.codeStart){return 1;}
		
		if(codeStop<other.codeStop){return -1;}
		if(codeStop>other.codeStop){return 1;}

		if(exonLength<other.exonLength){return -1;}
		if(exonLength>other.exonLength){return 1;}
		
		if(strand<other.strand){return -1;}
		if(strand>other.strand){return 1;}
		
		if(id<other.id){return -1;}
		if(id>other.id){return 1;}
		
		if(!symbol.equals(other.symbol)){return symbol.compareTo(other.symbol);}
		return mrnaAcc.compareTo(other.mrnaAcc);
	}
	
	public boolean isIdenticalTo(Gene other){
		if(chromosome!=other.chromosome){return false;}
		if(strand!=other.strand){return false;}
		if(txStart!=other.txStart){return false;}
		if(txStop!=other.txStop){return false;}
		if(codeStart!=other.codeStart){return false;}
		if(codeStop!=other.codeStop){return false;}
		if(exonLength!=other.exonLength){return false;}
//		if(pseudo!=other.pseudo || untranslated!=other.untranslated){return false;}
		if(exons==null){
			if(other.exons!=null){return false;}
		}else{
			if(other.exons==null || (other.exons.length!=exons.length)){return false;}
			for(int i=0; i<exons.length; i++){
				Exon e1=exons[i], e2=other.exons[i];
				if(e1.a!=e2.a || e1.b!=e2.b){return false;}
				//assert(e1.equals(e2));
				//if(!e1.equals(e2)){return false;}
			}
		}
		return true;
	}
	
	@Override
	public boolean equals(Object other){
		return equals((Gene)other);
	}
	
	public boolean equals(Gene other){//TODO check this
		return this==other ? true : compareTo(other)==0;
	}
	
	@Override
	public int hashCode(){
		int xor=txStart^(Integer.rotateLeft(txStop, 16));
		xor^=(chromosome<<20);
		xor^=strand;
		return xor;
	}
	
	public int translationStart(){ //TODO Make into a field
		return (exons==null || exons.length==0) ? codeStart : exons[0].a;
	}
	
	public int translationStop(){ //TODO Make into a field
		return (exons==null || exons.length==0) ? codeStop : exons[exons.length-1].b;
	}
	
	public int codeStartStrandCompensated(){ //TODO Make into a field
		return strand==Shared.PLUS ? codeStart : codeStop;
	}
	
	public int codeStopStrandCompensated(){ //TODO Make into a field
		return strand==Shared.PLUS ? codeStop : codeStart;
	}
	
	public int translationStartStrandCompensated(){ //TODO Make into a field
		if(strand==Shared.PLUS){
			return (exons==null || exons.length==0) ? codeStart : exons[0].a;
		}else{
			return (exons==null || exons.length==0) ? codeStop : exons[exons.length-1].b;
		}
	}
	
	public int translationStopStrandCompensated(){ //TODO Make into a field
		if(strand==Shared.PLUS){
			return (exons==null || exons.length==0) ? codeStop : exons[exons.length-1].b;
		}else{
			return (exons==null || exons.length==0) ? codeStart : exons[0].a;
		}
	}
	
	public int exonStartStrandCompensated(int exNum){
		if(strand==Shared.PLUS){
			return (exons==null || exons.length==0) ? codeStart : exons[exNum].a;
		}else{
			return (exons==null || exons.length==0) ? codeStop : exons[exons.length-exNum-1].b;
		}
	}
	
	public int exonStopStrandCompensated(int exNum){
		if(strand==Shared.PLUS){
			return (exons==null || exons.length==0) ? codeStop : exons[exNum].b;
		}else{
			return (exons==null || exons.length==0) ? codeStart : exons[exons.length-exNum-1].a;
		}
	}

	public int findClosestExon(int a, int b) {
		if(exons==null || exons.length==0){return 0;}
		int best=Integer.MAX_VALUE;
		int exnum=-1;
		for(int i=0; i<exons.length; i++){
			Exon e=exons[i];
			int x=calcDistance(a, b, e.a, e.b);
			if(x<best){
				best=x;
				if(strand==Shared.PLUS){exnum=i;}
				else{exnum=exons.length-i-1;}
			}
		}
		assert(exnum>=0);
		return exnum;
	}
	
	
	/** Calculates the minimal distance between two ranges: (a1, b1) and (a2, b2). */
	public static final int calcDistance(int a1, int b1, int a2, int b2){
		assert(a1<=b1);
		assert(a2<=b2);
		int r;
		if(b1>=a2 && b2>=a1){r=0;}
		else if(a1>b2){r=a1-b2;}
		else{r=a2-b1;}
		assert(r>=0) : r;
		return r;
	}
	

	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/** Transcription start position */
	public final int txStart;
	/** Transcription end position */
	public final int txStop;
	/** Coding region start */
	public final int codeStart;
	/** Coding region end */
	public final int codeStop;
	
	/** Length of transcribed area */
	public final int txLength;
	
	/** Length of coding area */
	public final int codeLength;
	
	/** Length of exons (summed) */
	public final int exonLength;
	
	/** Length of exonic coding region */
	public final int exonCodeLength;
	
	/** Number of amino acids (excluding stop codon) */
	public final int aaLength;
	
	public final int utrLength5prime;
	public final int utrLength3prime;

	/** Reference sequence chromosome or scaffold */
	public final byte chromosome;
	/** + or - for strand */
	public final byte strand;
	/** ? */
	public final byte cdsStartStat;
	/** ? */
	public final byte cdsEndStat;
	
	public final boolean readCorrectly;
	
	/** Array of exons used by this gene */
	public final Exon[] exons;
	
	/** Exon frame {0,1,2}, or -1 if no frame for exon */
	public final byte[] exonFrames;
	
	/** Name of gene (usually transcript_id from GTF) */
	public final String mrnaAcc;
	
	/** Protein accession */
	public final String proteinAcc;
	
	/** Alternate name (e.g. gene_id from GTF) */
	public final String symbol;
	
	public final String description;
	
	public final String fullDescription;
	
	public final byte primarySource;
	
	/* CCDS file format:
	 * chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type */
		
	/* CCDS format stuff */

//	public final String nc_accession;
	public final byte status;
	public final byte completeness;
	public final int id;
	

	public final boolean untranslated;
	public final boolean pseudo;
	public final boolean valid;
	
	public static final String[] sourceCodes={
		"seqGene", "knownGene", "refGene", "unionGene",
		"reserved1", "reserved2", "reserved3", "reserved4"
	};
	
	/** Index with cdsStartStat and cdsEndStat */
	public static final String[] endStatCodes={"none", "unk", "incmpl", "cmpl"};
	
	public static final String[] statusCodes={
		"Unknown","Reviewed","Validated","Provisional","Predicted","Inferred","Public"
		
//		"Public", "Reviewed, update pending", "Reviewed, withdrawal pending",
//		"Withdrawn", "Withdrawn, inconsistent annotation",
//		"Under review, withdrawal", "Under review, update",
		
	};
	
	public static final String[] completenessCodes={
		"Unknown","Complete5End","Complete3End","FullLength","IncompleteBothEnds","Incomplete5End","Incomplete3End","Partial"
	};
	
	
	/** Index with chromosome number */
	public static final String[] chromCodes={"A", "1", "2", "3", "4", "5", "6",
		"7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
		"19", "20", "21", "22", "X", "Y", "M", "U"};
	
	private static final int NEAR=Data.NEAR;

	public static final byte STAT_UNKNOWN=0;
	public static final byte STAT_REVIEWED=1;
	public static final byte STAT_VALIDATED=2;
	public static final byte STAT_PROVISIONAL=3;
	public static final byte STAT_PREDICTED=4;
	public static final byte STAT_INFERRED=5;
	public static final byte STAT_PUBLIC=6;

}
