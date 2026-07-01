package var2;

import java.io.PrintStream;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Generates neural-network training vectors from a label-annotated VCF (the output of
 * {@link LabelVCF}).  This is the depth-independent half of the pipeline: the NN target rides
 * with the variant (the INFO key NNL, computed once on the deepest data), while the feature
 * vector is computed from whatever depth THIS VCF represents.  The same variant therefore gets
 * the same target but a different vector at 10x, 40x, 160x.
 *
 * The VCF is streamed one line at a time: each variant is parsed with {@link VcfToVar#fromVCF},
 * its label read from NNL, and the {@link VectorUMP45} feature vector emitted followed by the
 * target.  Variants whose label falls in the deadzone (0.5 +/- deadzone, default the 0.4-0.6
 * band) are excluded so ambiguous examples never reach the net.  Streaming keeps memory flat
 * regardless of file size, and the vector math is reused verbatim from VectorUMP45 — only the
 * label SOURCE differs from VcfToTrainingVectors (read, not recomputed).
 *
 * Dataset-level statistics (pairing rate, quality/mapq/readlen averages, ploidy) are read from
 * the VCF ## header lines, the same ones VcfLoader consumes.
 *
 * Usage: java var2.VectorsFromLabeledVCF in=labeled.vcf.gz out=vectors.tsv [ref=ref.fa.gz]
 *        [deadzone=0.1] [ploidy=N] [platform=illumina] [includescore=f]
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class VectorsFromLabeledVCF {

	public static void main(String[] args){
		String inFile=null, refFile=null, outFile=null;
		double deadzone=0.1;        // half-width of the excluded band around 0.5 (0.1 => skip 0.4-0.6)
		int ploidyOverride=-1;      // if >0, overrides the ##ploidy from the header
		HashSet<String> allowedChroms=null; // if set, only these CHROM values are vectorized

		for(String arg : args){
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in")){inFile=b;}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("deadzone")){deadzone=Double.parseDouble(b);}
			else if(a.equals("ploidy")){ploidyOverride=Integer.parseInt(b);}
			else if(a.equals("includescore")){VectorUMP45.includeScore=Parse.parseBoolean(b);}
			else if(a.equals("scaffolds") || a.equals("chroms")){
				allowedChroms=new HashSet<String>();
				for(String c : b.split(",")){allowedChroms.add(c.trim());}
			}
			else if(a.equals("platform")){
				if(b.equalsIgnoreCase("illumina")){VectorUMP45.platform=VectorUMP45.PLATFORM_ILLUMINA;}
				else if(b.equalsIgnoreCase("pacbio")){VectorUMP45.platform=VectorUMP45.PLATFORM_PACBIO;}
				else if(b.equalsIgnoreCase("nanopore") || b.equalsIgnoreCase("ont")){VectorUMP45.platform=VectorUMP45.PLATFORM_NANOPORE;}
				else if(b.equalsIgnoreCase("roche") || b.equalsIgnoreCase("sbx")){VectorUMP45.platform=VectorUMP45.PLATFORM_ROCHE;}
				else{VectorUMP45.platform=Integer.parseInt(b);}
			}
		}

		assert(inFile!=null) : "Missing in= (labeled VCF)";
		assert(outFile!=null) : "Missing out= (vector TSV)";
		assert(deadzone>=0 && deadzone<0.5) : "deadzone must be in [0,0.5): "+deadzone;

		PrintStream out=System.err;

		// Scaffold map for VcfToVar.fromVCF / makeVector.  Per-variant features (HMP, CED, etc.)
		// come from the extended INFO already present in the VCF, so the header map suffices when
		// no reference is supplied.
		ScafMap scafMap;
		if(refFile!=null){
			out.println("Loading reference: "+refFile);
			scafMap=ScafMap.loadReference(refFile, true);
		}else{
			out.println("Loading scaffold map from VCF header: "+inFile);
			scafMap=ScafMap.loadVcfHeader(inFile);
		}

		// Dataset-level stats (overwritten from ## header lines as we stream).
		double pairingRate=0, totalQualityAvg=0, totalMapqAvg=0, readLengthAvg=0;
		int ploidy=2;

		final double dzLo=0.5-deadzone, dzHi=0.5+deadzone;

		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(inFile, FileFormat.TXT, null, true, false));
		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		bb.append("#dims\t").append(VectorUMP45.DIMS).append('\t').append(1).nl();
		bsw.print(bb);

		long written=0, excluded=0, positive=0, negative=0, nolabel=0, filtered=0;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]=='#'){
				// Capture dataset stats from ## header lines (key=value).
				String h=new String(line);
				int eq=h.indexOf('=');
				if(eq>0){
					String k=h.substring(0, eq), val=h.substring(eq+1).trim();
					if(k.equalsIgnoreCase("##properPairRate")){pairingRate=Double.parseDouble(val);}
					else if(k.equalsIgnoreCase("##totalQualityAvg")){totalQualityAvg=Double.parseDouble(val);}
					else if(k.equalsIgnoreCase("##mapqAvg")){totalMapqAvg=Double.parseDouble(val);}
					else if(k.equalsIgnoreCase("##readLengthAvg")){readLengthAvg=Double.parseDouble(val);}
					else if(k.equalsIgnoreCase("##ploidy")){ploidy=Integer.parseInt(val);}
				}
				continue;
			}

			// Read the NNL label from the INFO column.
			String s=new String(line);
			if(allowedChroms!=null){
				int t=s.indexOf('\t');
				String chrom=(t>0 ? s.substring(0, t) : s);
				if(!allowedChroms.contains(chrom)){filtered++; continue;}
			}
			int idx=s.indexOf(";NNL=");
			if(idx<0){nolabel++; continue;}
			int lstart=idx+5;
			int lend=s.indexOf(';', lstart);
			if(lend<0){lend=s.length();}
			double label=Double.parseDouble(s.substring(lstart, lend));

			// Deadzone exclusion (ambiguous labels never train).
			if(label>dzLo && label<dzHi){excluded++; continue;}

			// Parse the variant and compute its feature vector.
			Var v;
			try{
				v=VcfToVar.fromVCF(line, scafMap, true, true);
			}catch(Exception e){
				nolabel++; continue;
			}
			if(v==null){nolabel++; continue;}
			v.calcCoverage(scafMap);

			final int useploidy=(ploidyOverride>0 ? ploidyOverride : ploidy);
			float[] vec=VectorUMP45.makeVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, useploidy, scafMap);

			bb.clear();
			for(int i=0; i<vec.length; i++){
				if(i>0){bb.tab();}
				bb.append(vec[i], 6);
			}
			bb.tab().append(label, 6).nl();
			bsw.print(bb);

			written++;
			if(label>=0.5){positive++;}else{negative++;}
		}
		bf.close();
		//TODO: Possible bug [var2/VectorsFromLabeledVCF#001] (LOW, dropped-wire) - poisonAndWait()'s errorState
		//return is discarded; a write failure -> truncated vector TSV while the tool exits 0. Capture + exit(1).
		bsw.poisonAndWait();

		out.println("Vectors written: "+written+" (positive: "+positive+", negative: "+negative+")");
		out.println("Excluded (deadzone "+dzLo+"-"+dzHi+"): "+excluded);
		out.println("Filtered (scaffold not in include list): "+filtered);
		out.println("Skipped (no label / parse fail): "+nolabel);
		out.println("Output: "+outFile);
	}
}
