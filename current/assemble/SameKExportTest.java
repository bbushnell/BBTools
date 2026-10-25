package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import map.ObjectIntMap;
import parse.LineParser1;

/** Independent small-fixture read counting and full exported-vector verification.
 * @author Fischl */
final class SameKExportTest {

	/** Reads FASTA fixtures only; no native read-count table is reused by this oracle. */
	static void verify(final String reads, final String directory, final int k){
		final ObjectIntMap<String> counts=new ObjectIntMap<String>(String.class);
		for(Contig c : TadpoleGraph.readContigs(reads)){
			final String read=text(c.bases);
			for(int p=0; p+k<=read.length(); p++){
				final String key=canonical(read.substring(p, p+k));
				//Keep this literal count oracle independent of ObjectIntMap.increment.
				final int old=counts.get(key);
				counts.put(key, old<0 ? 1 : old+1);
			}
		}
		long observations=0, squares=0;
		final String[] keys=counts.keys();
		final int[] values=counts.values();
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=null){observations+=values[i]; squares+=(long)values[i]*values[i];}
		}
		final ArrayList<Contig> contigs=TadpoleGraph.readContigs(directory+"/raw.fa");
		final ByteFile features=ByteFile.makeByteFile(directory+"/features.tsv", true);
		final ByteFile candidates=ByteFile.makeByteFile(directory+"/candidates.tsv", true);
		final ByteFile controls=ByteFile.makeByteFile(directory+"/controls.tsv", true);
		final LineParser1 fp=new LineParser1('\t'), cp=new LineParser1('\t');
		int candidateCount=0, controlCount=0;
		boolean failed=false;
		try{
			final String[] metadata=text(features.nextLine()).split(" ");
			check(metadata.length==4 && metadata[0].equals("#schema="+SameKJoinFeatures.VERSION), "Missing feature schema.");
			check(Long.parseLong(metadata[1].substring(7))==counts.size(), "Native unique depth population differs from scalar reads.");
			check(Double.parseDouble(metadata[2].substring(13))==observations, "Native observation sum differs from scalar reads.");
			check(Double.parseDouble(metadata[3].substring(8))==squares, "Native squared-depth sum differs from scalar reads.");
			features.nextLine();
			check(text(candidates.nextLine()).equals(SameKJoinDiagnostic.HEADER), "Changed candidate schema.");
			check(text(controls.nextLine()).equals(SameKJoinExport.CONTROL_HEADER), "Wrong virtual-control schema.");
			final float[] expected=new float[27];
			for(byte[] row=features.nextLine(); row!=null; row=features.nextLine()){
				if(row[0]=='#'){continue;}
				fp.set(row);
				check(fp.terms()==29, "Wrong feature width.");
				final boolean control=fp.termEquals("control", 0);
				check(control || fp.termEquals("candidate", 0), "Unknown feature population.");
				final byte[] context=(control ? controls : candidates).nextLine();
				check(context!=null, "Missing matching context row.");
				cp.set(context);
				check(cp.parseInt(0)==fp.parseInt(1) && cp.parseInt(1)==k, "Mismatched feature/context ID or K.");
				final String a, b;
				if(control){
					controlCount++;
					final String whole=text(contigs.get(cp.parseInt(2)).bases);
					final int cut=cp.parseInt(3);
					a=whole.substring(0, cut); b=whole.substring(cut-k+1);
					check(cp.parseInt(4)==a.length() && cp.parseInt(5)==b.length(), "Wrong virtual arm lengths.");
					check(text(cp.parseByteArray(6)).equals(a.substring(a.length()-k)) &&
							text(cp.parseByteArray(7)).equals(b.substring(0, k)), "Wrong control endpoint words.");
					check(text(cp.parseByteArray(8)).equals(a.substring(Math.max(0, a.length()-k+1-200))) &&
							text(cp.parseByteArray(9)).equals(b.substring(0, Math.min(b.length(), k-1+200))), "Wrong bounded control contexts.");
				}else{
					candidateCount++;
					final String source=text(contigs.get(cp.parseInt(2)).bases), dest=text(contigs.get(cp.parseInt(4)).bases);
					a=cp.parseInt(3)==0 ? rc(source) : source;
					b=cp.parseInt(5)==1 ? rc(dest) : dest;
				}
				final long asum=sum(a, k, counts), bsum=sum(b, k, counts);
				final int an=a.length()-k+1, bn=b.length()-k+1;
				final double am=control ? asum/(double)an : asum/(float)an;
				final double bm=control ? bsum/(double)bn : bsum/(float)bn;
				SameKJoinFeatures.fill(expected, k, neighborhood(a.substring(a.length()-k), counts),
						neighborhood(b.substring(0, k), counts), "ACGT".indexOf(a.charAt(a.length()-k-1)),
						"ACGT".indexOf(b.charAt(k)), am, bm, an, bn,
						observations/(double)counts.size(), squares/(double)observations);
				for(int i=0; i<27; i++){
					check(fp.parseFloat(i+2)==expected[i], "Scalar feature disagreement: kind="+(control ? "control" : "candidate")+
							" event="+fp.parseInt(1)+" field="+i);
				}
			}
			check(candidates.nextLine()==null && controls.nextLine()==null, "Context rows lack matching feature vectors.");
			check(candidateCount>0 && controlCount==candidateCount/9, "Fixture expected default one control per nine candidates.");
		}finally{failed=features.close(); failed|=candidates.close(); failed|=controls.close();}
		check(!failed, "Fixture input failed to close.");
		System.err.println("SAME_K_EXPORT_TEST_PASS k="+k+" unique="+counts.size()+" observations="+observations+
				" squares="+squares+" candidates="+candidateCount+" controls="+controlCount);
	}

	/** Scalar substring scan, independent of the rolling Kmer control sampler. */
	private static long sum(final String sequence, final int k, final ObjectIntMap<String> counts){
		long sum=0;
		for(int p=0; p+k<=sequence.length(); p++){
			final int count=counts.get(canonical(sequence.substring(p, p+k)));
			check(count>0, "Raw formed contig has no scalar read evidence.");
			sum+=count;
		}
		return sum;
	}

	/** Constructs all oriented neighbors with strings, not the production packed-Kmer helpers. */
	private static int[] neighborhood(final String word, final ObjectIntMap<String> counts){
		final int[] result=new int[9]; result[0]=Math.max(0, counts.get(canonical(word)));
		for(int b=0; b<4; b++){
			result[1+b]=Math.max(0, counts.get(canonical(word.substring(1)+"ACGT".charAt(b))));
			result[5+b]=Math.max(0, counts.get(canonical("ACGT".charAt(b)+word.substring(0, word.length()-1))));
		}
		return result;
	}

	private static String canonical(final String word){final String reverse=rc(word); return word.compareTo(reverse)<=0 ? word : reverse;}
	private static String text(final byte[] bases){check(bases!=null, "Unexpected end of fixture input."); return new String(bases, StandardCharsets.US_ASCII);}
	private static String rc(final String word){
		final StringBuilder result=new StringBuilder(word.length());
		for(int i=word.length()-1; i>=0; i--){result.append("TGCA".charAt("ACGT".indexOf(word.charAt(i))));}
		return result.toString();
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
