package var2;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import parse.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Represents a single line from a VCF (Variant Call Format) file.
 * Provides parsing, manipulation, and splitting capabilities for genetic variants.
 * Supports VCF standard fields including scaffold, position, reference, alternate,
 * quality, filters, and sample data. Handles complex variants by splitting into
 * simpler components and canonical representations.
 *
 * @author Brian Bushnell
 */
public class VCFLine implements Comparable<VCFLine>, Cloneable {
	
	/**
	 * Constructs a VCFLine from individual field components.
	 *
	 * @param scaf_ Scaffold/chromosome name
	 * @param pos_ Position on the scaffold (1-based)
	 * @param id_ Variant identifier
	 * @param ref_ Reference allele sequence
	 * @param alt_ Alternate allele sequence
	 * @param qual_ Quality score
	 * @param filter_ Filter status
	 * @param info_ INFO field data
	 * @param format_ FORMAT field specification
	 * @param type_ Variant type constant from Var class
	 * @param samples_ Sample data list
	 */
	public VCFLine(String scaf_, int pos_, byte[] id_, byte[] ref_, byte[] alt_, double qual_, 
			byte[] filter_, byte[] info_, byte[] format_, int type_, ArrayList<byte[]> samples_) {
		scaf=scaf_;
		pos=pos_;
		id=id_;
		ref=ref_;
		reflen=ref[0]=='.' ? 0 : ref.length;
		alt=alt_;
		qual=qual_;
		filter=filter_;
		info=info_;
		format=format_;
		type=type_;
		if(samples_!=null) {
			samples.addAll(samples_);
		}
		hashcode=hash();
	}
	
	/**
	 * Constructs a VCFLine by parsing a tab-delimited VCF line.
	 * Automatically determines variant type and processes sample data.
	 * Applies canonical trimming and caching if enabled.
	 * @param line Raw VCF line as byte array
	 */
	public VCFLine(byte[] line) {
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		scaf=new String(line, a, b-a, StandardCharsets.US_ASCII);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		pos=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		id=line[a]=='.' ? DOT : Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		if(b<=a+1){ref=Var.AL_MAP[line[a]];}
		else{ref=Arrays.copyOfRange(line, a, b);}
		
		reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		if(b<=a+1){alt=Var.AL_MAP[line[a]];}
		else{alt=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		qual=(line[a]=='.' && b==a+1 ? 40 : Parse.parseDouble(line, a, b));
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		filter=Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		info=Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		int TYP=Tools.indexOfDelimited(info, "TYP=", 0, (byte)';');
		type=(TYP>=0 ? Var.typeInitialArray[info[TYP+4]] : type_old());
		assert(type>=0) : type+", "+TYP+"\n"+new String(info);
		
		while(b<line.length && line[b]!='\t'){b++;}
//		assert(b>a) : "Missing field 8: "+new String(line);
		if(b>a){//LoFreq does not produce Info field
			format=Arrays.copyOfRange(line, a, b);
			b++;
			a=b;
		}
		
		while(b<line.length){
			while(b<line.length && line[b]!='\t'){b++;}
			if(b<=a){
				break;
			}
			byte[] sample=Arrays.copyOfRange(line, a, b);
			samples.add(sample);
			b++;
			a=b;
		}
		
		if(TRIM_TO_CANONICAL){trimToCanonical();}
//		assert(false) : "a="+a+", b="+b+", "+(b<=a+1)+", ref="+new String(ref)+"\n"+new String(line)+"\n"+this;
		
		hashcode=hash();
		if(AUTOCACHE){cache();}
	}
	
	/** Converts this VCFLine to a Var object for internal processing.
	 * @return Var representation of this variant */
	public Var toVar(){
		return makeVar(info, alt);
	}
	
	/**
	 * Creates a Var object from VCF INFO field and alternate allele.
	 * Parses structured INFO data including coverage, quality metrics,
	 * strand bias information, and statistical measures.
	 *
	 * @param info INFO field containing semicolon-delimited key=value pairs
	 * @param alt Alternate allele sequence
	 * @return Var object with parsed statistics and metadata
	 */
	public static Var makeVar(byte[] info, byte[] alt){
		int a=0, b=0;
		
		//SN=0;STA=547693;STO=547694;TYP=SUB;
		assert(Tools.startsWith(info, "SN", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 0: "+new String(info);
		int scaf=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 1: "+new String(info);
		int start=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 2: "+new String(info);
		int stop=Parse.parseInt(info, a, b);
		b++;
		a=b;

//		assert(Tools.startsWith(info, "TYP", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 3: "+new String(info);
		final int type=Var.typeInitialArray[info[a]];
		assert(type>=0) : type+new String(info);
//		if(Tools.contains(info, SUB, a)){type=Var.SUB;}
//		else if(Tools.contains(info, DEL, a)){type=Var.DEL;}
//		else if(Tools.contains(info, INS, a)){type=Var.INS;}
//		else if(Tools.contains(info, NOCALL, a)){type=Var.NOCALL;}
//		else{assert(false) : new String(info);}
		b++;
		a=b;
		
		//R1P=20;R1M=29;R2P=25;R2M=19;
		assert(Tools.startsWith(info, "R1P", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 4: "+new String(info);
		int r1p=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 5: "+new String(info);
		int r1m=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 6: "+new String(info);
		int r2p=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 7: "+new String(info);
		int r2m=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		//AD=2;DP=24;MCOV=0;PPC=0;
		assert(Tools.startsWith(info, "AD=", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 8: "+new String(info);
//		int ad=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 9: "+new String(info);
		int cov=Parse.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 10: "+new String(info);
		int mcov=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "PPC", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 11: "+new String(info);
		int pc=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		//AF=0.0833;RAF=0.0833;LS=280;
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 12: "+new String(info);
//		double af=Parse.parseDouble(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "RAF", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 13: "+new String(info);
		double raf=Parse.parseDouble(info, a, b);
		b++;
		a=b;
		
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 14: "+new String(info);
		long ls=Parse.parseLong(info, a, b);
		b++;
		a=b;
		
		//MQS=86;MQM=43;BQS=64;BQM=32;
		assert(Tools.startsWith(info, "MQS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 15: "+new String(info);
		long mqs=Parse.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 16: "+new String(info);
		int mqm=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "BQS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 17: "+new String(info);
		long bqs=Parse.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 18: "+new String(info);
		int bqm=Parse.parseInt(info, a, b);
		b++;
		a=b;

		//EDS=18;EDM=9;IDS=1984;IDM=992;NVC=0;FLG=0;
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 19: "+new String(info);
		long eds=Parse.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 20: "+new String(info);
		int edm=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "IDS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 21: "+new String(info);
		long ids=Parse.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 22: "+new String(info);
		int idm=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "NVC", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 23: "+new String(info);
		int nvc=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "FLG", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 24: "+new String(info);
		int flg=Parse.parseInt(info, a, b);
		b++;
		a=b;
		
//		//CED=13;HMP=1;SB=0.9980;DP=24;DP4=22,0,2,0
//		assert(Tools.startsWith(info, "CED", a));
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 23: "+new String(info);
////		int ced=Parse.parseInt(info, a, b);
//		b++;
//		a=b;
//		
//		//HMP=2;AF=0.989;
//		assert(Tools.startsWith(info, "HMP", a));
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 24: "+new String(info);
////		int hmp=Parse.parseInt(info, a, b);
//		b++;
//		a=b;
//		
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 25: "+new String(info);
////		double sb=Parse.parseDouble(info, a, b);
//		b++;
//		a=b;
//
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 26: "+new String(info);
////		String dp4=new String(info, a, b-a);
//		b++;
//		a=b;
		
		if(type==Var.DEL || type==Var.INS){
			if(alt.length<=1){alt=Var.AL_0;}
			else if(alt.length==2){alt=Var.AL_MAP[alt[1]];}
			else{alt=Arrays.copyOfRange(alt, 1, alt.length);}
		}
		
		//GT:DP:AD:AF	1:94:93:0.989

		Var v=new Var(scaf, start, stop, alt, type);
		v.r1plus=r1p;
		v.r1minus=r1m;
		v.r2plus=r2p;
		v.r2minus=r2m;
		v.properPairCount=pc;
		v.lengthSum=ls;
		v.mapQSum=mqs;
		v.mapQMax=mqm;
		v.baseQSum=bqs;
		v.baseQMax=bqm;
		v.endDistSum=eds;
		v.endDistMax=edm;
		v.idSum=ids;
		v.idMax=idm;
		v.nearbyVarCount=nvc;
		v.setFlagged(flg>0);
		v.revisedAlleleFraction=raf;
		v.setCoverage(cov, mcov);
//		v.homopolymerCount=hmp; //derived
		
		return v;
	}
	
	/**
	 * Splits complex variants into simpler components based on specified criteria.
	 * Can separate multi-allelic variants, complex variants, and multi-base substitutions.
	 *
	 * @param splitAlleles Split multi-allelic variants (comma-separated alternatives)
	 * @param splitComplex Split complex variants into substitution and indel components
	 * @param splitSubs Split multi-base substitutions into individual SNPs
	 * @return List of simplified VCFLine objects, or null if no splitting needed
	 */
	public ArrayList<VCFLine> split(boolean splitAlleles, boolean splitComplex, boolean splitSubs){
		assert(splitAlleles || splitComplex || splitSubs);
		if(isSimple()){return null;} //Should be true 90% of the time
		if(isJunction()) {assert(false) : "TODO";}
		
		if(!isMulti()){splitAlleles=false;}
//		if(type!=Var.SUB || alt.length<2){splitSubs=false;}
		
		assert(!splitComplex) : "splitComplex function is written, but not integrated.";
		
		if(splitAlleles){
			ArrayList<VCFLine> list0=splitAlleles();
			if(!splitSubs || list0==null){
				if(SORT && list0!=null){Collections.sort(list0);}
				if(CONDENSE){condense(list0);}
				return list0;
			}
			ArrayList<VCFLine> list2=new ArrayList<VCFLine>();
			for(VCFLine line0 : list0){
				ArrayList<VCFLine> list1=null;
				if(line0.type==Var.SUB && line0.alt.length>1){
					list1=line0.splitSubs();
					assert(list1!=null) : "\n"+this+"\n"+line0+"\n";
				}
				if(list1!=null){list2.addAll(list1);}
				else{list2.add(line0);}
			}
			if(SORT && list2!=null){Collections.sort(list2);}
			if(CONDENSE){condense(list2);}
			return list2;
		}
		
		if(type!=Var.SUB || alt.length<2){splitSubs=false;}
		if(!splitSubs || isMulti()){return null;}//Can't split subs with commas
		
		ArrayList<VCFLine> list=splitSubs();
		return list;
	}
	
	/** Split into one line per allele */
	private ArrayList<VCFLine> splitAlleles(){
		assert(isMulti()) : this;
		if(Tools.indexOf(alt, ',')<0){
			assert(false);
			return null;
		}
		String alleles=new String(alt);
		String[] split=alleles.split(",");
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(split.length);
//		System.err.println(this);
//		System.err.println(new String(ref)+", "+new String(alt)+", "+type());
//		System.err.println(this);
		
		final String[] splitInfo=(info==null || !SPLIT_INFO ? null : new String(info).split(";"));
		
		for(int i=0; i<split.length; i++){
			final String s=split[i];
			VCFLine line=this.clone();
			line.alt=s.getBytes();

//			System.err.println(line);
//			System.err.println(line.type());
			line.recalc();
//			System.err.println(line);
//			System.err.println(line.type());
			if(!line.isRef()){
				if(TRIM_TO_CANONICAL){line.trimToCanonical();}
				if(SPLIT_INFO){line.info=splitInfo(splitInfo, i, split.length);}
//				System.err.println(line);
				list.add(line);
			}
		}
		return list.size()==0 ? null : list;
	}
	
	/** Splits multi-base substitutions into SNPs.
	 * Discards resultant ref "SNPs" */ 
	private ArrayList<VCFLine> splitSubs(){
		assert(type==Var.SUB);
		assert(alt.length>1);
		assert(Tools.indexOf(alt, ',')<0) : toString();
		assert(alt.length==ref.length) : this;
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(alt.length);
		for(int i=0; i<alt.length; i++){
//			VCFLine line=this.clone();
//			line.alt=Var.AL_MAP[alt[i]];
//			line.ref=Var.AL_MAP[ref[i]];
//			line.pos+=i;
//			line.reflen=1;
//			line.rehash();
			VCFLine line=new VCFLine(scaf, pos+i, id, Var.AL_MAP[ref[i]], Var.AL_MAP[alt[i]], qual, filter, info, format, type, samples);
			if(!line.isRef() && (Var.CALL_NOCALL || !line.isNocall())){
				if(TRIM_TO_CANONICAL){line.trimToCanonical();}
				list.add(line);
			}
		}
		return list.size()==0 ? null : list;
	}
	
	/** Splits non-length-neutral lines into sub+del or sub+ins.
	 * Discards resultant ref "SUBs"
	 * Since alignment information is no longer present, this will do a bad job usually,
	 * unless alignment is performed. 
	 * It's OK for easy things like a 1bp substitution plus a single deletion or insertion. */
	private ArrayList<VCFLine> splitComplex(){
		assert(type==Var.COMPLEX);
		assert(alt.length>1);
		assert(Tools.indexOf(alt, ',')<0) : toString();
		assert(alt.length!=ref.length && alt.length>1 && ref.length>1) : this;
		
		if(TRIM_TO_CANONICAL){this.trimToCanonical();}
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(2);

		final int prefixLen=Tools.min(alt.length, ref.length)-1;
		final byte[] prefixA=prefix(alt, prefixLen);
		final byte[] prefixR=prefix(ref, prefixLen);
		final byte[] suffixA=suffix(alt, alt.length-prefixLen);
		final byte[] suffixR=suffix(ref, ref.length-prefixLen);
		VCFLine sub=new VCFLine(scaf, pos, id, prefixR, prefixA, qual, filter, info, format, type, samples);
		VCFLine indel=new VCFLine(scaf, pos, id, suffixR, suffixA, qual, filter, info, format, type, samples);
		
		assert(sub.isSub());
		if(alt.length<ref.length){//Sub + del
			assert(indel.isDel());
		}else if(alt.length>ref.length){//Sub + ins
			assert(indel.isIns());
		}else{
			assert(false) : this;
			throw new RuntimeException("Unreachable");
		}
		
		if(!sub.isRef() && (Var.CALL_NOCALL || !sub.isNocall())){list.add(sub);}
		list.add(indel);
		return list;
	}
	
	/** Attempts to split comma-delimited info fields by allele
	 * This will not always be correct as some info fields are supposed to contain commas */
	private byte[] splitInfo(String[] splitInfo, int alleleNum, int alleles){
		if(splitInfo==null){return null;}
		assert(alleles>1);
		assert(alleleNum>=0 && alleleNum<alleles);
		ByteBuilder bb=new ByteBuilder();
		for(String part : splitInfo){
			if(bb.length()>0){bb.append(';');}
			if(part.indexOf(',')<0){
				bb.append(part);
			}else{
				String[] splitEquals=part.split("=");
				String[] splitComma=splitEquals[1].split(",");
				if(splitComma.length==alleles){
					bb.append(splitEquals[0]).append('=').append(splitComma[alleleNum]);
				}else{
					bb.append(part);
				}
			}
		}
		return bb.toBytes();
	}
	
	/** Assumes list is sorted; removes duplicates */
	private static int condense(ArrayList<VCFLine> list){
		if(list==null || list.size()<2){return 0;}
		int removed=0;
		VCFLine prev=list.get(0);
		for(int i=1; i<list.size(); i++){
			VCFLine line=list.get(i);
			if(prev.equals(line)){
				list.set(i, null);
				removed++;
			}else{
				prev=line;
			}
		}
		if(removed>0){Tools.condenseStrict(list);}
		return removed;
	}
	
	@Override
	public VCFLine clone(){
		try {
			return (VCFLine)super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	
	/** Trim matching affixes for a canonical representation */
	private int trimToCanonical(){
		if(isMulti() || alt.length<2 || ref.length<2 /*|| ref.length==alt.length*/){return 0;}
		assert(Tools.indexOf(alt, ',')<0) : this;
		
		int s=trimSuffix();
		int p=trimPrefix();
		int delta=s+p;
		if(delta>0){recalc();}
		assert(!isRef()) : this;
		return delta;
	}
	
	/**
	 * Removes matching prefix bases from reference and alternate alleles.
	 * Updates position to reflect trimmed sequence start.
	 * @return Number of prefix bases removed
	 */
	private int trimPrefix(){
		int prefix=0;
		for(int rpos=0, apos=0; rpos<ref.length-1 && apos<alt.length-1 && ref[rpos]==alt[apos]; rpos++, apos++){
			prefix++;
		}
		if(prefix>0){
			pos+=prefix;
			reflen=ref.length-prefix;
			int altlen=alt.length-prefix;
			ref=(reflen==1 ? Var.AL_MAP[ref[ref.length-1]] : Arrays.copyOfRange(ref, prefix, ref.length));
			alt=(altlen==1 ? Var.AL_MAP[alt[alt.length-1]] : Arrays.copyOfRange(alt, prefix, alt.length));
			assert(alt.length>0  && ref.length>0) : this;
		}
		return prefix;
	}
	
	/** Removes matching suffix bases from reference and alternate alleles.
	 * @return Number of suffix bases removed */
	private int trimSuffix(){
		int suffix=0;
		for(int rpos=ref.length-1, apos=alt.length-1; rpos>0 && apos>0 && ref[rpos]==alt[apos]; rpos--, apos--){
			suffix++;
		}
		if(suffix>0){
			reflen=ref.length-suffix;
			int altlen=alt.length-suffix;
			ref=(reflen==1 ? Var.AL_MAP[ref[0]] : Arrays.copyOf(ref, reflen));
			alt=(altlen==1 ? Var.AL_MAP[alt[0]] : Arrays.copyOf(alt, altlen));
			assert(alt.length>0  && ref.length>0) : this;
		}
		return suffix;
	}
	
	/**
	 * Extracts a prefix of specified length from a byte array.
	 * @param array Source array
	 * @param len Length of prefix to extract
	 * @return Prefix bytes, optimized for single-byte case
	 */
	private byte[] prefix(byte[] array, int len){
		assert(len>0);
		assert(array.length>len);
		byte[] prefix=(len==1 ? Var.AL_MAP[array[0]] : len==array.length ? array : Arrays.copyOf(array, len));
		return prefix;
	}
	
	/**
	 * Extracts a suffix of specified length from a byte array.
	 * @param array Source array
	 * @param len Length of suffix to extract
	 * @return Suffix bytes, optimized for single-byte case
	 */
	private byte[] suffix(byte[] array, int len){
		assert(len>0);
		assert(array.length>len);
		byte[] suffix=(len==1 ? Var.AL_MAP[array[array.length-1]] : len==array.length ? array : Arrays.copyOfRange(array, array.length-len, array.length));
		return suffix;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b){
		return equals((VCFLine)b);
	}
	
	/**
	 * Tests equality with another VCFLine using hash codes and comparison.
	 * @param b VCFLine to compare against
	 * @return true if lines represent identical variants
	 */
	public boolean equals(VCFLine b){
		return hashcode==b.hashcode && compareTo(b)==0;
	}

	/** Computes hash code based on scaffold, position, length, and alternate allele.
	 * @return Hash code for this variant */
	private int hash(){
		return scaf.hashCode()^Integer.rotateLeft(pos, 9)^Integer.rotateRight(pos+ref.length, 9)^Var.hash(alt);
	}

//	private void rehash(){
//		hashcode=hash();
//	}
	
	@Override
	public int hashCode(){
		return hashcode;
	}
	
	/** Generates a unique key for indexing based on position and hash.
	 * @return 62-bit key combining position, hash, and scaffold information */
	public long toKey() {
		long key=Long.rotateLeft(pos, 31)^Long.rotateRight(hashcode, 10)^scaf.hashCode();
		return key&0x3FFFFFFFFFFFFFFFL;
	}
	
	@Override
	public int compareTo(VCFLine v){
		ScafMap map=ScafMap.defaultScafMap();
		assert(map!=null);
		int scafnum1=map.getScaffold(scaf).number;
		int scafnum2=map.getScaffold(v.scaf).number;
		if(scafnum1!=scafnum2){return scafnum1-scafnum2;}
		if(pos!=v.pos){return pos-v.pos;}
		final int typeA=type(), typeB=v.type();
		if(typeA!=typeB){return typeA-typeB;}
		int stop1=pos+reflen(), stop2=v.pos+reflen();
		if(stop1!=stop2){return stop1-stop2;}
		return compare(alt, v.alt);
	}
	
	/**
	 * Compares two byte arrays lexicographically with length precedence.
	 * @param a First array
	 * @param b Second array
	 * @return Comparison result favoring shorter arrays when different lengths
	 */
	public int compare(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a.length!=b.length){return b.length-a.length;}
		for(int i=0; i<a.length; i++){
			byte ca=a[i], cb=b[i];
			if(ca!=cb){return ca-cb;}
		}
		return 0;
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		return toText(bb).toString();
	}
	
	/**
	 * Appends VCF-formatted representation to a ByteBuilder.
	 * @param bb ByteBuilder to append to
	 * @return The modified ByteBuilder
	 */
	public ByteBuilder toText(ByteBuilder bb){
		bb.append(scaf).append('\t');
		bb.append(pos).append('\t');
		bb.append(id).append('\t');
		bb.append(ref).append('\t');
		bb.append(alt).append('\t');
		bb.append(qual, 2).append('\t');
		bb.append(filter).append('\t');
		bb.append(info).append('\t');
		if(format!=null){
			bb.append(format);
		}
		for(byte[] sample : samples){
			bb.tab().append(sample);
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the length of the reference allele */
	public int reflen(){return ref.length;}
	/** Gets the length of the alternate allele */
	public int readlen(){return alt.length;}
	
	/**
	 * Determines variant type using legacy classification logic.
	 * Categorizes based on allele lengths and content analysis.
	 * @return Variant type constant from Var class
	 */
	public int type_old(){
		if(alt!=null && Tools.indexOf(alt, ',')>=0){return Var.MULTI;}
		final int reflen=reflen(), readlen=readlen();
		if(readlen!=reflen && readlen!=0 && reflen!=0){return Var.COMPLEX;}
		if(reflen<readlen){return Var.INS;}
		if(reflen>readlen){return Var.DEL;}
		for(byte b : alt){
			if(b!='N'){return Var.SUB;}
		}
		return Var.NOCALL;
	}
	
	/** Gets the variant type constant */
	public int type(){return type;}
	
	/** Tests if alternate allele matches reference (no variant) */
	public boolean isRef(){
		return Tools.equals(alt, ref);
	}
	
	/** Tests if variant represents a structural variant junction */
	public boolean isJunction(){
		return type==Var.LJUNCT || type==Var.RJUNCT || type==Var.BJUNCT;
	}
	
	/** Tests if variant is an insertion or deletion */
	public boolean isIndel(){
		return type==Var.INS || type==Var.DEL;
	}
	
	/** Tests if variant is a substitution */
	public boolean isSub(){
		return type==Var.SUB;
	}
	
	/** Tests if variant is a deletion */
	public boolean isDel(){
		return type==Var.DEL;
	}
	
	/** Tests if variant is an insertion */
	public boolean isIns(){
		return type==Var.INS;
	}
	
	/** Tests if variant represents a no-call (N bases) */
	public boolean isNocall(){
		return type==Var.NOCALL;
	}
	
	/** Tests if variant has multiple alternate alleles */
	public boolean isMulti(){
		return type==Var.MULTI;
	}
	
	/** Tests if variant is a complex rearrangement */
	public boolean isComplex(){
		return type==Var.COMPLEX;
	}
	
	/** Tests if variant is a simple type (SUB, DEL, INS, or NOCALL) */
	public boolean isSimple(){
		//return type<=Var.DEL;
		return type==Var.SUB || type==Var.DEL || type==Var.INS || type==Var.NOCALL;
	}
	
	/** Applies string caching to reduce memory usage for repeated sequences */
	void cache(){
//		assert(false) : AUTOCACHE;
		id=cache(id);
		if(ref.length<5){ref=cache(ref);}
		if(alt.length<5){alt=cache(alt);}
		filter=cache(filter);
		ref=cache(ref);
		format=cache(format);
	}
	
	/**
	 * Caches byte arrays to reduce memory usage for repeated sequences.
	 * @param line Byte array to cache
	 * @return Cached instance or original if cache limit exceeded
	 */
	static byte[] cache(byte[] line){
		if(line==null){return line;}
		String s=new String(line);
		byte[] old=cache.get(s);
		if(old!=null){return old;}
		if(cache.size()>20000){return line;}
		synchronized(cache){
			old=cache.get(s);
			if(old!=null){return old;}
			cache.put(s, line);
			return line;
		}
	}
	
	/**
	 * Converts string to cached byte array.
	 * @param s String to cache
	 * @return Cached byte array representation
	 */
	static byte[] cache(String s){
		if(s==null){return null;}
		byte[] old=cache.get(s);
		if(old!=null){return old;}
		synchronized(cache){
			old=cache.get(s);
			if(old!=null){return old;}
			old=s.getBytes();
			cache.put(s, old);
			return old;
		}
	}
	
	/** 0-based */
	public int start(){
		return pos-1;
	}
	
	/** 0-based */
	public int stop(){
		return Tools.max(start(), pos+reflen-2);
	}
	
	/** Recalculates variant type and hash code after sequence modification */
	private void recalc(){
		type=type_old();
		hashcode=hash();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Scaffold/chromosome name */
	public final String scaf;
	/** Position on scaffold (1-based) */
	public int pos;
	/** Variant identifier */
	public byte[] id;
	/** Reference allele sequence */
	public byte[] ref;
	/** Length of reference allele */
	public int reflen;
	/** Alternate allele sequence */
	public byte[] alt;
	/** Quality score */
	public double qual;
	/** Filter status */
	public byte[] filter;
	/** INFO field containing variant annotations */
	public byte[] info;
	/** FORMAT field specifying sample data format */
	public byte[] format;
	/** Precomputed hash code for this variant */
	public int hashcode;
	/** Variant type constant from Var class */
	public int type;
	/** Sample data columns */
	public ArrayList<byte[]> samples=new ArrayList<byte[]>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String-to-byte array cache for memory optimization */
	public static HashMap<String, byte[]> cache=new HashMap<String, byte[]>(99997);
	
	/** Whether to automatically apply caching during construction */
	static boolean AUTOCACHE=false;
	/** Whether to trim variants to canonical representation */
	static boolean TRIM_TO_CANONICAL=true;
	/** Whether to sort split variant lists */
	static boolean SORT=true;
	/** Whether to remove duplicates from split variant lists */
	static boolean CONDENSE=true;
	/** Whether to split INFO fields by allele for multi-allelic variants */
	static boolean SPLIT_INFO=true;
	
	private static final byte[] NOCALL=cache("NOCALL");
	private static final byte[] SUB=cache("SUB");
	private static final byte[] DEL=cache("DEL");
	private static final byte[] INS=cache("INS");
	private static final byte[] LJUNCT=cache("LJUNCT");
	private static final byte[] RJUNCT=cache("RJUNCT");
	private static final byte[] BJUNCT=cache("BJUNCT");
	private static final byte[] MULTI=cache("MULTI");
	private static final byte[] COMPLEX=cache("COMPLEX");
	private static final byte[] DOT=cache(".");
	private static final byte[] PASS=cache("PASS");
	private static final byte[] FAIL=cache("FAIL");
	private static final byte[] FORMAT=cache("GT:DP:AD:AF:SC:PF");
	
	static{
		cache(Var.AL_0);
		cache(Var.AL_A);
		cache(Var.AL_C);
		cache(Var.AL_G);
		cache(Var.AL_T);
		cache(Var.AL_N);
	}

}
