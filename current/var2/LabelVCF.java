package var2;

import java.io.PrintStream;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.TextFile;
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
 * GIAB variant may legitimately target >1.0.
 *
 * The deadzone / blacklist decision (a label too close to 0.5) is deliberately left to the
 * downstream vector generator, so the deadzone width stays tunable without re-running this
 * program; here we only count how many labels fall in 0.4-0.6 for diagnostics.
 *
 * GIAB membership and high-confidence membership are computed with the SAME canonical key
 * logic that {@link VcfToTrainingVectors} uses (POS/REF/ALT with indel leading-base stripping),
 * so the GIAB VCF and the candidate VCF must be normalized the same way (e.g. both run through
 * filtervcf normalize splitalleles).  The GIAB VCF passed here must NOT be bed-restricted —
 * the BED is supplied separately and membership in it is what separates cat 1/3 from cat 2/4.
 *
 * Usage: java var2.LabelVCF in=allvariants.vcf.gz giab=giab_norm.vcf.gz bed=highconf.bed \
 *        out=labeled.vcf.gz [ref=ref.fa.gz] [divisor=55] [phredcap=60]
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class LabelVCF {

	public static void main(String[] args){
		String inFile=null, giabFile=null, bedFile=null, refFile=null, outFile=null;
		double giabHC=0.85, giabLC=0.75;     // GIAB offsets (in-bed, off-bed)
		double divisor=55;                   // non-GIAB phred divisor
		double hcPenalty=0.10, lcPenalty=0.05; // non-GIAB penalties (in-bed, off-bed)
		double phredCap=60;                  // clamp phred before applying formulas

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
		}

		assert(inFile!=null) : "Missing in= (candidate VCF)";
		assert(giabFile!=null) : "Missing giab= (GIAB truth VCF; normalized, NOT bed-restricted)";
		assert(bedFile!=null) : "Missing bed= (GIAB high-confidence BED)";
		assert(outFile!=null) : "Missing out= (labeled VCF)";
		assert(divisor>0) : "divisor must be positive: "+divisor;
		assert(phredCap>0) : "phredcap must be positive: "+phredCap;
		assert(giabLC<=giabHC) : "giablc ("+giabLC+") should not exceed giabhc ("+giabHC+")";

		PrintStream out=System.err;

		// Scaffold map: from the reference if given, else from the candidate VCF header (lighter,
		// avoids loading the whole reference just to map contig names to numbers).
		ScafMap scafMap;
		if(refFile!=null){
			out.println("Loading reference: "+refFile);
			scafMap=ScafMap.loadReference(refFile, true);
		}else{
			out.println("Loading scaffold map from VCF header: "+inFile);
			scafMap=ScafMap.loadVcfHeader(inFile);
		}

		// GIAB truth keys (NOT bed-restricted): membership separates cat 1/2 from cat 3/4.
		out.println("Loading GIAB truth keys: "+giabFile);
		HashSet<String> giab=loadKeys(giabFile, scafMap);
		out.println("GIAB variant keys: "+giab.size());

		// High-confidence BED: in-bed => cat 1 or 3, off-bed => cat 2 or 4.
		out.println("Loading high-confidence BED: "+bedFile);
		BedMask bed=new BedMask(bedFile);
		out.println("BED intervals: "+bed.intervalsLoaded()+" across "+bed.scaffolds()+" scaffolds");

		// Stream candidate VCF, splice NNL/NNC into INFO, emit.
		out.println("Labeling candidates: "+inFile);
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(inFile, FileFormat.TXT, null, true, false));
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		long total=0, cat1=0, cat2=0, cat3=0, cat4=0, blacklist=0, unkScaf=0;
		boolean headerInjected=false;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]=='#'){
				// Inject NNL/NNC INFO definitions just before the #CHROM column header.
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

			// Data line: split into VCF columns. CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE...]
			String[] f=new String(line).split("\t");
			if(f.length<8){ // malformed; pass through unchanged
				bb.clear(); bb.append(line).nl(); bsw.print(bb);
				continue;
			}
			String chrom=f[0];
			int pos=Integer.parseInt(f[1]);
			String ref=f[3];
			String alt=f[4];
			double phred=(f[5].length()==1 && f[5].charAt(0)=='.') ? 0 : Double.parseDouble(f[5]);
			if(phred>phredCap){phred=phredCap;}
			if(phred<0){phred=0;}

			String key=makeKeyFromVcfFields(chrom, pos, ref, alt, scafMap);
			if(key==null){unkScaf++;}
			boolean inGiab=(key!=null && giab.contains(key));
			boolean inBed=bed.contains(chrom, pos);

			double label;
			int cat;
			if(inGiab && inBed){label=giabHC+phred/100.0; cat=1; cat1++;}
			else if(inGiab){label=giabLC+phred/100.0; cat=2; cat2++;}
			else if(inBed){label=Math.max(0.0, phred/divisor-hcPenalty); cat=3; cat3++;}
			else{label=Math.max(0.0, phred/divisor-lcPenalty); cat=4; cat4++;}

			if(label>0.4 && label<0.6){blacklist++;}

			// Rebuild the line, appending ;NNL=...;NNC=... to the INFO column (index 7).
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
	 * Loads a VCF as a set of canonical variant keys (CHROM/POS/REF/ALT with indel leading-base
	 * stripping).  Identical key construction to {@link VcfToTrainingVectors} so the GIAB key set
	 * matches candidate keys built from the same scaffold map.  Works with any VCF format.
	 */
	private static HashSet<String> loadKeys(String fname, ScafMap scafMap){
		HashSet<String> set=new HashSet<String>();
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("#")){continue;}
			String[] fields=line.split("\t", 6);
			if(fields.length<5){continue;}

			String chrom=fields[0];
			int pos=Integer.parseInt(fields[1]);
			String ref=fields[3];
			String alt=fields[4];

			// Handle multi-allelic (comma-separated ALT)
			String[] alts=alt.split(",");
			for(String a : alts){
				String key=makeKeyFromVcfFields(chrom, pos, ref, a.trim(), scafMap);
				if(key!=null){set.add(key);}
			}
		}
		tf.close();
		return set;
	}

	/**
	 * Creates a canonical variant key from VCF fields for matching, converting VCF coordinates
	 * to the internal 0-based representation used by VcfToVar.  Returns null if the scaffold is
	 * not in the map (e.g. an off-target contig absent from GIAB's reference).
	 */
	private static String makeKeyFromVcfFields(String chrom, int pos, String ref, String alt, ScafMap scafMap){
		int scafNum=scafMap.getNumber(chrom);
		if(scafNum<0){return null;}

		int refLen=ref.length();
		int altLen=alt.length();

		int start, stop;
		String normalizedAlt;

		if(refLen!=altLen && altLen>0 && refLen>0){
			// Indel: strip leading base
			normalizedAlt=alt.substring(1);
			start=pos; // stays 1-based here, converted via refLen below
			refLen--;
			if(refLen==0 && normalizedAlt.isEmpty()){
				return null; // degenerate single-base insertion that got stripped
			}
		}else{
			normalizedAlt=alt;
			start=pos-1; // convert to 0-based
		}
		stop=start+refLen;

		return scafNum+":"+start+":"+stop+":"+normalizedAlt.toUpperCase();
	}
}
