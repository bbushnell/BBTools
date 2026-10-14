package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;
import ukmer.Kmer;

/** Scalar fixtures for the shared feature schema and internal-path control sampler.
 * @author Fischl */
final class SameKJoinFeatureTest {

	/** Invoked by the existing samekjointest launcher on both compiled classpaths. */
	static int run(){
		checks=0;
		normalization();
		serialization();
		sampling();
		for(int k : new int[]{3, 31, 32, 33, 62, 64, 95, 96, 124}){controls(k);}
		return checks;
	}

	/** Native decimal parsing must recover the exact floats consumed by inference. */
	private static void serialization(){
		final Random random=new Random(821);
		final ByteBuilder row=new ByteBuilder();
		final LineParser1 parser=new LineParser1('\t');
		final float[] vector=new float[27];
		for(int trial=0; trial<100; trial++){
			for(int i=0; i<vector.length; i++){
				vector[i]=SameKJoinFeatures.depth(random.nextInt(Integer.MAX_VALUE));
			}
			vector[0]=0; vector[1]=0.125f;
			row.clear().append("candidate\t1");
			SameKJoinExport.appendVector(row, vector);
			parser.set(row.toBytes());
			check(parser.terms()==29, "Feature export lost a field.");
			for(int i=0; i<27; i++){check(parser.parseFloat(i+2)==vector[i], "Float feature changed in text roundtrip.");}
		}
		boolean failed=false;
		vector[0]=Float.MIN_VALUE;
		try{SameKJoinExport.appendVector(row.clear(), vector);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "An unrepresentable future feature was silently rounded to zero.");
	}

	/** Opposite traversal swaps endpoints and complements/reverses both neighborhoods. */
	private static void normalization(){
		final int[] a={20, 0, 19, 2, 1, 17, 3, 0, 1}, b={18, 1, 0, 3, 14, 0, 2, 20, 1};
		final int[] before=a.clone();
		final float[] forward=new float[27], reverse=new float[27];
		SameKJoinFeatures.fill(forward, 96, a, b, 0, 3, 17.25, 19.5, 105, 300, 12, 45);
		SameKJoinFeatures.fill(reverse, 96, reverse(b), reverse(a), 0, 3, 19.5, 17.25, 300, 105, 12, 45);
		check(Arrays.equals(forward, reverse), "Partner/strand reversal changed the model inputs.");
		check(Arrays.equals(a, before), "Feature sorting mutated raw nucleotide-indexed counts.");
		check(SameKJoinFeatures.depth(0)==0 && SameKJoinFeatures.depth(1)==0.125f, "Absent and singleton depth were conflated.");
		check(SameKJoinFeatures.depth(Integer.MAX_VALUE)>3, "High depths were clipped.");
		for(int i=0; i<27; i++){check(Float.isFinite(forward[i]) && forward[i]>=0, "Invalid normalized field "+i);}
		final float[] expected=new float[27];
		// a sorts after b because current depth20>18. Spell all fields independently.
		final double[] raw={18, 20, 2, 1, 0, 14, 3, 1, 0, 14, 19.5, 300,
				20, 19, 2, 1, 0, 17, 3, 1, 0, 17, 17.25, 105, 96, 12, 45};
		for(int i=0; i<27; i++){
			expected[i]=(float)((i==11 || i==23 ? 0.0625*Math.log(raw[i]) :
					i==24 ? 0.125*Math.log(raw[i]) : 0.125*Math.log1p(raw[i]))/Math.log(2));
			check(Math.abs(forward[i]-expected[i])<1e-7, "Feature position/order disagreement at "+i);
		}
		final SameKJoinFeatures.Moments moments=new SameKJoinFeatures.Moments();
		check(moments.uniqueMean()==0 && moments.observationMean()==0, "Empty table moments are not a valid no-op.");
		moments.add(1); moments.add(3);
		check(moments.unique==2 && moments.uniqueMean()==2 && moments.observationMean()==2.5,
				"Unique-weighted and observation-weighted means were confused.");
		moments.add(Integer.MAX_VALUE);
		check(moments.observationMean()>2e9, "Count square overflowed before conversion to double.");
		for(int mode=0; mode<7; mode++){
			boolean failed=false;
			try{
				if(mode<3){SameKJoinFeatures.depth(mode==0 ? -1 : mode==1 ? Double.NaN : Double.POSITIVE_INFINITY);}
				else if(mode==3){moments.add(0);}
				else{
					SameKJoinFeatures.fill(new float[mode==4 ? 26 : 27], 96, a, b,
							mode==5 ? -1 : 0, 3, 1, 1, mode==6 ? 0 : 1, 1, 1, 1);
				}
			}catch(IllegalArgumentException expectedFailure){failed=true;}
			check(failed, "Invalid feature input accepted, mode="+mode);
		}
	}

	/** Full, empty, huge and seeded subset draws preserve range, size and uniqueness. */
	private static void sampling(){
		check(SameKAcceptedControls.sample(0, 0, 7).length==0, "Empty population sampling failed.");
		final long[] full=SameKAcceptedControls.sample(15, 15, 7);
		for(int i=0; i<full.length; i++){check(full[i]==i, "Complete population lost a cut.");}
		for(long population : new long[]{15, 1000, Long.MAX_VALUE}){
			final long[] a=SameKAcceptedControls.sample(population, 10, 59);
			check(Arrays.equals(a, SameKAcceptedControls.sample(population, 10, 59)), "Seeded sampling is not reproducible.");
			for(int i=0; i<a.length; i++){check(a[i]>=0 && a[i]<population && (i==0 || a[i]>a[i-1]), "Invalid sampled cut.");}
		}
	}

	/** Brute-force substring counts independently check every virtual-arm mean and edge. */
	private static void controls(final int k){
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		final Random random=new Random(1731+k);
		for(int length : new int[]{k+4, k+5, k+27}){
			final byte[] bases=new byte[length];
			for(int i=0; i<length; i++){bases[i]=AminoAcid.numberToBase[random.nextInt(4)];}
			contigs.add(new Contig(bases, contigs.size()));
		}
		final ObjectIntMap<String> table=new ObjectIntMap<String>(String.class);
		for(Contig c : contigs){
			for(int p=0; p+k<=c.length(); p++){
				table.put(canonical(new String(c.bases, p, k, StandardCharsets.US_ASCII)), 1+random.nextInt(90));
			}
		}
		final TadpoleGraph.Counts evidence=new TadpoleGraph.Counts(){
			@Override public int count(final Kmer word){return table.get(canonical(word.toString()));}
		};
		final long[] seen={0};
		final SameKAcceptedControls.Result result=new SameKAcceptedControls(contigs, evidence, k).collect(1000, 7,
				new SameKAcceptedControls.Sink(){
			@Override
			public void accept(int id, int cut, Contig c, double am, double bm, int[] ac, int[] bc, int ai, int bi){
				seen[0]++;
				check(cut>=k+2 && cut<=c.length()-3 && id>0, "Near-K virtual control escaped the formation guard.");
				long left=0, right=0;
				for(int p=0; p+k<=c.length(); p++){
					final int value=table.get(canonical(new String(c.bases, p, k, StandardCharsets.US_ASCII)));
					if(p+k<=cut){left+=value;}else{right+=value;}
				}
				check(am==left/(double)(cut-k+1) && bm==right/(double)(c.length()-cut), "Virtual-arm depth means disagree with scalar scan.");
				check(ai==AminoAcid.baseToNumber[c.bases[cut-k-1]] && bi==AminoAcid.baseToNumber[c.bases[cut+1]],
						"Actual inward continuation has an off-by-one or complement error.");
				check(ac[0]==table.get(canonical(new String(c.bases, cut-k, k, StandardCharsets.US_ASCII))) &&
						bc[0]==table.get(canonical(new String(c.bases, cut-k+1, k, StandardCharsets.US_ASCII))), "Endpoint Kmer depth mismatch.");
				final float[] f=new float[27], r=new float[27];
				SameKJoinFeatures.fill(f, k, ac, bc, ai, bi, am, bm, cut-k+1, c.length()-cut, 20, 40);
				SameKJoinFeatures.fill(r, k, reverse(bc), reverse(ac), 3-bi, 3-ai, bm, am, c.length()-cut, cut-k+1, 20, 40);
				check(Arrays.equals(f, r), "Control features are not reversal invariant at K="+k);
			}
		});
		check(result.population==24 && result.sampled==24 && seen[0]==24 && result.probability()==1,
				"Eligible virtual split population should be 0+1+23.");
		final SameKAcceptedControls.Result empty=new SameKAcceptedControls(contigs, evidence, k).collect(0, 7,
				new SameKAcceptedControls.Sink(){
			@Override public void accept(int i, int p, Contig c, double a, double b, int[] x, int[] y, int ai, int bi){
				throw new AssertionError("Zero requested controls emitted a cut.");
			}
		});
		check(empty.sampled==0 && empty.population==24 && empty.probability()==0, "Zero sampling lost population metadata.");
	}

	/** Scalar reverse-complement transformation of the unsorted nine-count schema. */
	private static int[] reverse(final int[] counts){
		final int[] result=new int[9]; result[0]=counts[0];
		for(int i=0; i<4; i++){result[1+i]=counts[8-i]; result[5+i]=counts[4-i];}
		return result;
	}

	/** Test-only exact string key, independent of Kmer's canonical packed representation. */
	private static String canonical(final String word){
		final StringBuilder reverse=new StringBuilder(word.length());
		for(int i=word.length()-1; i>=0; i--){reverse.append("TGCA".charAt("ACGT".indexOf(word.charAt(i))));}
		final String rc=reverse.toString();
		return word.compareTo(rc)<=0 ? word : rc;
	}

	private static void check(final boolean condition, final String message){checks++; if(!condition){throw new AssertionError(message);}}
	private static int checks;
}
