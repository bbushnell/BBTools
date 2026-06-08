package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Annotates a CallVariants VCF with a neural-network training-target label (INFO key NNL),
 * producing a depth-independent "labeled master" truth set.
 *
 * The label is computed from two axes:
 *   1) GIAB truth membership  (in the GIAB VCF or not), and
 *   2) GIAB high-confidence region membership (inside the high-confidence BED or not),
 * combined with the candidate's own CallVariants QUAL score (which is read directly from
 * the VCF QUAL column — it is the same value CallVariants writes as the FORMAT SC field).
 *
 * Four label categories (written as INFO key NNC):
 *   cat 1  GIAB,  in-bed   : label = giabHC + phred/100          (default 0.85 + phred/100)
 *   cat 2  GIAB,  off-bed  : label = giabLC + phred/100          (default 0.75 + phred/100)
 *   cat 3  !GIAB, in-bed   : label = max(0, phred/divisor - hcPenalty)  (default phred/55 - 0.10)
 *   cat 4  !GIAB, off-bed  : label = max(0, phred/divisor - lcPenalty)  (default phred/55 - 0.05)
 *
 * The phred score is clamped to phredCap (default 60) before the formula is applied; this
 * matches the empirical CallVariants ceiling (~Q60) and bounds the non-GIAB tail.  Targets
 * are intentionally NOT clamped to 1.0 — the NN output function is unbounded, so a strong
 * GIAB variant may legitimately target >1.0.  When phredCap=0, all phred modulation is
 * suppressed and the labels become binary (GIAB offset vs 0).
 *
 * The deadzone / blacklist decision (a label too close to 0.5) is deliberately left to the
 * downstream vector generator, so the deadzone width stays tunable without re-running this
 * program; here we only count how many labels fall in 0.4-0.6 for diagnostics.
 *
 * GIAB membership uses VCFLine equals/hashCode — the same identity contract as CompareVCF —
 * so concordance with comparevcf intersection is exact.  When ref= is given, indels are
 * left-aligned for normalization; this matches the split=t normalize pipeline.
 *
 * Usage: java var2.LabelVCF in=allvariants.vcf.gz giab=giab_norm.vcf.gz bed=highconf.bed \
 *        out=labeled.vcf.gz [ref=ref.fa.gz] [divisor=55] [phredcap=60] [normalize=t]
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class LabelVCF {

	public static void main(String[] args){
		String inFile=null, giabFile=null, bedFile=null, refFile=null, outFile=null;
		double giabHC=0.85, giabLC=0.75;
		double divisor=55;
		double hcPenalty=0.10, lcPenalty=0.05;
		double phredCap=60;
		boolean normalize=false;

		for(String arg : args){
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in")){inFile=b;}
			else if(a.equals("giab") || a.equals("truth")){giabFile=b;}
			else if(a.equals("bed")){bedFile=b;}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("giabhc")){giabHC=Double.parseDouble(b);}
			else if(a.equals("giablc")){giabLC=Double.parseDouble(b);}
			else if(a.equals("divisor")){divisor=Double.parseDouble(b);}
			else if(a.equals("hcpenalty")){hcPenalty=Double.parseDouble(b);}
			else if(a.equals("lcpenalty")){lcPenalty=Double.parseDouble(b);}
			else if(a.equals("phredcap")){phredCap=Double.parseDouble(b);}
			else if(a.equals("normalize") || a.equals("leftalign") || a.equals("norm")){normalize=Parse.parseBoolean(b);}
			else if(a.equals("split") || a.equals("splitalleles")){/* always split truth */}
		}

		assert(inFile!=null) : "Missing in= (candidate VCF)";
		assert(giabFile!=null) : "Missing giab= (GIAB truth VCF; normalized, NOT bed-restricted)";
		assert(bedFile!=null) : "Missing bed= (GIAB high-confidence BED)";
		assert(outFile!=null) : "Missing out= (labeled VCF)";
		assert(divisor>0) : "divisor must be positive: "+divisor;
		assert(phredCap>=0) : "phredcap must be non-negative: "+phredCap;
		assert(giabLC<=giabHC) : "giablc ("+giabLC+") should not exceed giabhc ("+giabHC+")";

		if(normalize && refFile==null){
			throw new RuntimeException("normalize requires ref=");
		}

		PrintStream out=System.err;

		if(refFile!=null){
			out.println("Loading reference: "+refFile);
			ScafMap.loadReference(refFile, null, null, true);
		}
		ScafMap scafMap;
		if(ScafMap.defaultScafMap()!=null){
			scafMap=ScafMap.defaultScafMap();
		}else{
			out.println("Loading scaffold map from VCF header: "+inFile);
			scafMap=ScafMap.loadVcfHeader(inFile);
			ScafMap.setDefaultScafMap(scafMap, inFile);
		}

		out.println("Loading GIAB truth: "+giabFile);
		HashSet<VCFLine> giab=loadTruthSet(giabFile, scafMap, normalize);
		out.println("GIAB variants: "+giab.size());

		out.println("Loading high-confidence BED: "+bedFile);
		BedMask bed=new BedMask(bedFile);
		out.println("BED intervals: "+bed.intervalsLoaded()+" across "+bed.scaffolds()+" scaffolds");

		out.println("Labeling candidates: "+inFile);
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(inFile, FileFormat.TXT, null, true, false));
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		long total=0, cat1=0, cat2=0, cat3=0, cat4=0, blacklist=0, unkScaf=0;
		boolean headerInjected=false;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]=='#'){
				if(!headerInjected && Tools.startsWith(line, "#CHROM")){
					bb.clear();
					bb.append("##INFO=<ID=NNL,Number=1,Type=Float,Description=\"UMP45 NN training target label\">").nl();
					bb.append("##INFO=<ID=NNC,Number=1,Type=Integer,Description=\"UMP45 label category: 1=GIAB-HC,2=GIAB-LC,3=nonGIAB-HCregion,4=nonGIAB-LCregion\">").nl();
					bsw.print(bb);
					headerInjected=true;
				}
				bb.clear();
				bb.append(line).nl();
				bsw.print(bb);
				continue;
			}

			VCFLine vline;
			try{
				vline=new VCFLine(line);
			}catch(Throwable e){
				bb.clear(); bb.append(line).nl(); bsw.print(bb);
				continue;
			}

			if(scafMap.getScaffold(vline.scaf)==null){unkScaf++;}
			if(normalize){vline.leftAlign(refBasesFor(vline.scaf, scafMap));}

			boolean inGiab=giab.contains(vline);
			boolean inBed=bed.contains(vline.scaf, vline.pos);

			double phred=vline.qual;
			if(phred>phredCap){phred=phredCap;}
			if(phred<0){phred=0;}

			double label;
			int cat;
			if(inGiab && inBed){label=giabHC+phred/100.0; cat=1; cat1++;}
			else if(inGiab){label=giabLC+phred/100.0; cat=2; cat2++;}
			else if(inBed){label=Math.max(0.0, phred/divisor-hcPenalty); cat=3; cat3++;}
			else{label=Math.max(0.0, phred/divisor-lcPenalty); cat=4; cat4++;}

			if(label>0.4 && label<0.6){blacklist++;}

			String[] f=new String(line).split("\t");
			bb.clear();
			for(int i=0; i<f.length; i++){
				if(i>0){bb.tab();}
				bb.append(f[i]);
				if(i==7){bb.append(";NNL=").append(label, 5).append(";NNC=").append(cat);}
			}
			bb.nl();
			bsw.print(bb);
			total++;
		}
		bf.close();
		bsw.poisonAndWait();

		out.println("Labeled variants: "+total);
		out.println("  cat1 GIAB,  in-bed:   "+cat1);
		out.println("  cat2 GIAB,  off-bed:  "+cat2);
		out.println("  cat3 !GIAB, in-bed:   "+cat3);
		out.println("  cat4 !GIAB, off-bed:  "+cat4);
		out.println("  in deadzone 0.4-0.6:  "+blacklist);
		out.println("  unknown-scaffold:     "+unkScaf);
		out.println("Output: "+outFile);
	}

	/**
	 * Loads a truth VCF into a HashSet of VCFLine objects, using the same identity
	 * contract (equals/hashCode) as CompareVCF. Always splits multi-allelic variants.
	 * Optionally left-aligns indels for normalization.
	 */
	private static HashSet<VCFLine> loadTruthSet(String fname, ScafMap scafMap,
			boolean normalize){
		HashSet<VCFLine> set=new HashSet<VCFLine>();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TXT, null, true, false);
		ByteFile bf=ByteFile.makeByteFile(ff);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			VCFLine vline;
			try{
				vline=new VCFLine(line);
			}catch(Throwable e){continue;}
			if(scafMap.getScaffold(vline.scaf)==null){continue;}
			ArrayList<VCFLine> variants=vline.split(true, false, true);
			if(variants==null){
				if(normalize){vline.leftAlign(refBasesFor(vline.scaf, scafMap));}
				set.add(vline);
			}else{
				for(VCFLine v : variants){
					if(normalize){v.leftAlign(refBasesFor(v.scaf, scafMap));}
					set.add(v);
				}
			}
		}
		bf.close();
		return set;
	}

	private static byte[] refBasesFor(String scaf, ScafMap scafMap){
		ScafMap dm=ScafMap.defaultScafMap();
		ScafMap map=(dm!=null ? dm : scafMap);
		Scaffold sc=map.getScaffold(scaf);
		return sc==null ? null : sc.bases;
	}
}
