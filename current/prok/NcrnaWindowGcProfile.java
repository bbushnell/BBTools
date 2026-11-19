package prok;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Immutable, model-bound window-GC moments from the complete saved baseline.
 * Names, not row positions, bind the grader's output to the caller's library.
 * Known ACGT bases define GC, matching CaptureNcrnaWindows' saved observations.
 * Unknown window/model GC retains the unadjusted score rather than inventing a
 * rejection. A zero-SD profile admits exactly its observed GC. All fitting was
 * done upstream; no input window can change these moments.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWindowGcProfile{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaWindowGcProfile(final int models, final double weight_, final double cutoff_, final boolean range_){
		require(models>0 && Double.isFinite(weight_) && weight_>=0
			&& (range_ || Double.isFinite(cutoff_) && cutoff_>0),
			"GC scoring requires actual models, finite nonnegative weight and a positive SD cutoff or observed-range mode");
		weight=weight_;
		cutoff=cutoff_;
		range=range_;
		known=new boolean[models];
		mean=new double[models];
		sd=new double[models];
		minimum=new double[models];
		maximum=new double[models];
	}

	/** Cold metadata parsing deliberately uses round-trip Java doubles, including
	 * scientific notation. No strings or objects are allocated during scoring. */
	public static NcrnaWindowGcProfile load(final String path, final String[] names,
			final double weight, final double cutoff, final boolean range){
		require(path!=null && names!=null && names.length>0, "GC profiles must bind the actual nonempty caller library");
		final NcrnaWindowGcProfile result=new NcrnaWindowGcProfile(names.length, weight, cutoff, range);
		final ObjectIntMap<String> ordinals=new ObjectIntMap<String>(names.length, String.class);
		for(int model=0; model<names.length; model++){
			require(names[model]!=null, "Every GC profile needs a named actual caller model");
			final String name=RrnaResourceIO.first(names[model]);
			require(name.startsWith("cobalamin_RF00174_") && !ordinals.contains(name), "GC profile model identities must be unique RF00174 names");
			ordinals.put(name, model);
		}
		final LineParser1 p=new LineParser1('\t');
		final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
		final boolean[] seen=new boolean[names.length];
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==8, "Use the complete global baseline model-GC table");
				// The baseline writer uses BuildNcrnaReference.field, including
				// percent escaping. Decode once before binding a cold model row.
				p.setBounds(0);
				final int model=ordinals.get(decoder.decode(line, p.a(), p.b()));
				require(model>=0 && !seen[model], "A GC resource cannot add foreign models or repeat a caller model");
				seen[model]=true;
				count++;
				final long members=p.parseLong(1), missing=p.parseLong(2);
				require(members>=0 && missing>=0, "Known and unknown model-GC member counts cannot be negative");
				if(members==0){
					require(p.termEquals("NO_KNOWN_GC_MEMBERS", 7), "An unfitted model must remain explicitly unknown");
					for(int column=3; column<=6; column++){require(p.termEquals("NA", column), "Missing GC moments cannot silently become zero");}
					continue;
				}
				final double mean=Double.parseDouble(p.parseString(3)), sd=Double.parseDouble(p.parseString(4));
				final double min=Double.parseDouble(p.parseString(5)), max=Double.parseDouble(p.parseString(6));
				require(Double.isFinite(mean) && Double.isFinite(sd) && Double.isFinite(min) && Double.isFinite(max)
					&& min>=0 && min<=mean && mean<=max && max<=1 && sd>=0,
					"GC moments must be finite physical fractions with mean inside the observed range");
				require(sd==0 ? min==max && mean==min && p.termEquals("ZERO_SD", 7)
					: members>1 && max>min && p.termEquals("FIT", 7),
					"Zero-SD and fitted states must agree with the saved counts and actual range");
				result.known[model]=true;
				result.mean[model]=mean;
				result.sd[model]=sd;
				result.minimum[model]=min;
				result.maximum[model]=max;
				result.knownModels++;
			}
		}
		require(count==names.length, "Every actual caller model needs an explicit fitted or unknown GC row");
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Scoring           ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces every eligibility flag; adjusts only eligible, known-GC scores.
	 * Hard exclusions use a separate mask, never an infinity/NaN sentinel or a
	 * positive-score floor. A negative adjusted score is still a valid proposal.
	 * @return Number of models excluded by this window's GC gate */
	public int apply(final byte[] bases, final float[] scores, final boolean[] eligible){
		require(bases!=null && scores!=null && eligible!=null
			&& scores.length==known.length && eligible.length==known.length,
			"GC scoring must replace complete worker-owned score and eligibility vectors");
		int gc=0, acgt=0;
		for(final byte base:bases){
			switch(base|32){
				case 'c':
				case 'g':
					gc++;
					acgt++;
					break;
				case 'a':
				case 't':
					acgt++;
					break;
				default: break;
			}
		}
		final double fraction=acgt==0 ? Double.NaN : gc/(double)acgt;
		int rejected=0;
		for(int model=0; model<scores.length; model++){
			require(Float.isFinite(scores[model]) && scores[model]>=0, "GC penalties require the fresh nonnegative equal-bin score vector");
			eligible[model]=true;
			if(acgt==0 || !known[model]){continue;}
			final double difference=Math.abs(fraction-mean[model]);
			final double deviation=sd[model]==0 ? 0 : difference/sd[model];
			final boolean admitted=sd[model]==0 ? difference==0 : range
				? fraction>=minimum[model] && fraction<=maximum[model] : deviation<=cutoff;
			if(!admitted){
				eligible[model]=false;
				rejected++;
				continue;
			}
			final double adjusted=scores[model]-weight*deviation;
			require(Double.isFinite(adjusted) && adjusted>=-Float.MAX_VALUE && adjusted<=Float.MAX_VALUE,
				"A finite model-GC penalty must fit the caller's float score vector");
			scores[model]=(float)adjusted;
		}
		return rejected;
	}

	public int models(){return known.length;}
	public int knownModels(){return knownModels;}
	public double weight(){return weight;}
	public double cutoff(){return cutoff;}
	public boolean observedRange(){return range;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean[] known;
	private final double[] mean, sd, minimum, maximum;
	private final double weight, cutoff;
	private final boolean range;
	private int knownModels;
	public static final String HEADER="model\tknown_gc_members\tunknown_gc_members\tmean\tpopulation_sd\tminimum\tmaximum\tstate";
}
