package prok;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Immutable training-supported corrections on one loaded family's exact model axis.
 * Medians are signed caller errors in gene sense; centering subtracts them.
 * Workers borrow this resource and supply their own two-integer output scratch.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointModelOffsets{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads every model, including unsupported zero-offset rows, in library order. */
	public static NcrnaEndpointModelOffsets load(final String path, final NcrnaFamily family,
			final String member, final int minimumSupport){
		if(path==null || family==null || family.modelNames==null || family.library==null
				|| family.modelNames.length!=family.library.length || family.library.length==0
				|| member==null || minimumSupport<1){
			throw new IllegalArgumentException("Model offsets require a complete loaded family, member and positive support floor");
		}
		final NcrnaEndpointModelOffsets result=new NcrnaEndpointModelOffsets(family);
		final LineParser1 p=new LineParser1('\t');
		int model=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, NcrnaFixedEndpointTrials.OFFSET_HEADER, path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==6, "Model-offset rows require six fields: "+path);
				if(!p.termEquals(member, 0)){continue;}
				require(model<result.five.length && p.termEquals(NcrnaSupport.recordName(family.modelNames[model]), 1),
					"Offsets must retain the complete ordered model axis; model="+model+", file="+path);
				final long support=p.parseLong(2), five=p.parseLong(4), three=p.parseLong(5);
				require(support>=0 && p.parseInt(3)==minimumSupport
					&& five>=Integer.MIN_VALUE && five<=Integer.MAX_VALUE && three>=Integer.MIN_VALUE && three<=Integer.MAX_VALUE
					&& (support>=minimumSupport || five==0 && three==0),
					"Offset support/floor/range differs, or an unsupported model moves; model="+model+", file="+path);
				result.five[model]=(int)five;
				result.three[model]=(int)three;
				result.support[model]=support;
				model++;
			}
		}
		require(model==result.five.length, "Offsets cannot omit models; observed="+model+", expected="+result.five.length+", file="+path);
		return result;
	}

	private NcrnaEndpointModelOffsets(final NcrnaFamily family_){
		family=family_;
		five=new int[family.modelNames.length];
		three=new int[five.length];
		support=new long[five.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** A resource cannot be reinterpreted against another family's model order. */
	public boolean belongsTo(final NcrnaFamily other){return family==other;}
	public long trainingSupport(final int model){return support[model];}
	/** Signed training error on the declared sense-oriented end; centering subtracts it. */
	public int median(final int end, final int model){
		assert(end==0 || end==1) : "Model offsets have exactly the caller's 5prime/3prime axis; end="+end;
		return end==0 ? five[model] : three[model];
	}

	/** Writes the guarded center, or both original ends on invalid joint geometry.
	 * Physical fixed ends retain their original coordinates. A valid center is
	 * the fallback for unavailable features and rejected subsequent NN proposals.
	 */
	public boolean center(final int model, final int start, final int stop, final int fixedMask,
			final int contigLength, final int minimumLength, final int maximumLength, final int[] result){
		assert(model>=0 && model<five.length && result!=null && result.length==2 && (fixedMask&~3)==0
			&& contigLength>0 && minimumLength>0 && maximumLength>=minimumLength)
			: "The bound caller supplies its winning model, physical-end mask and joint length limits; model="+model;
		result[0]=start;
		result[1]=stop;
		final long nextStart=(long)start-((fixedMask&1)==0 ? five[model] : 0);
		final long nextStop=(long)stop-((fixedMask&2)==0 ? three[model] : 0);
		final long length=nextStop-nextStart+1;
		if(nextStart<0 || nextStop>=contigLength || length<minimumLength || length>maximumLength){return false;}
		result[0]=(int)nextStart;
		result[1]=(int)nextStop;
		return true;
	}

	private static void require(final boolean condition, final String why){
		if(!condition){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaFamily family;
	private final int[] five, three;
	private final long[] support;
}
