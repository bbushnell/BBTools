package prot;

import java.util.HashSet;
import java.util.Locale;

import prot.FamilyMatchInformation.View;
import structures.ByteBuilder;

/** Immutable stop-tail bonus formula for replay sweeps over one measured assignment.
 * The formula returns an absolute caller score contribution; callers must not multiply it
 * by the legacy tailbonus scale.
 * @author Yelan
 */
public final class FamilyBonusFormula{

	private FamilyBonusFormula(Form form_,Signal signal_,LengthMode lengthMode_,
			double a_,double b_,double x_,double y_,double z_,double scale_,double cap_,double outer_,
			double shortA_,double shortB_,double shortScale_,double shortOuter_,
			double i0_,double longI0_,double shortI0_,boolean automaticI0_,int boundMatrixScale_){
		form=form_;signal=signal_;lengthMode=lengthMode_;a=a_;b=b_;x=x_;y=y_;z=z_;scale=scale_;cap=cap_;i0=i0_;
		outer=outer_;shortA=shortA_;shortB=shortB_;shortScale=shortScale_;shortOuter=shortOuter_;
		longI0=longI0_;shortI0=shortI0_;
		automaticI0=automaticI0_;boundMatrixScale=boundMatrixScale_;
	}

	public boolean requiresInformation(){return signal!=null && signal.information;}
	public boolean requiresView(){return form==Form.SQRT_RATIO && longI0>0;}
	public boolean requiresFamilySetView(){return requiresView() || lengthMode==LengthMode.VIEW_STEP;}
	public boolean usesRawAlignmentScore(){return signal==Signal.BLOSUM_RAW;}

	/** Resolve an omitted raw divisor without mutating the requested formula.
	 * Rebinding checks every loaded matrix; an automatic scalar cannot represent
	 * different scales, even when an explicit-divisor override was requested.
	 */
	public FamilyBonusFormula bindRawScoreI0(final int matrixScale,final View view,
			final String binding,final boolean explicitOverride){
		if(!usesRawAlignmentScore()){return this;}
		require(matrixScale>0,"Bound matrix scale must be positive for "+binding);
		if(automaticI0){
			require(boundMatrixScale==0 || boundMatrixScale==matrixScale,
				"Automatic raw i0 requires the same matrix scale in every binding: "+binding);
			if(boundMatrixScale==0){
				return new FamilyBonusFormula(form,signal,lengthMode,a,b,x,y,z,scale,cap,outer,
					shortA,shortB,shortScale,shortOuter,
					matrixScale,longI0,shortI0,true,matrixScale);
			}
		}
		requireRawScoreI0(matrixScale,view,binding,explicitOverride);
		return this;
	}

	/** Raw alignment scores are emitted in the bound matrix's integer scale. */
	public void requireRawScoreI0(final int matrixScale,final View view,
			final String binding,final boolean explicitOverride){
		if(!usesRawAlignmentScore()){return;}
		require(matrixScale>0,"Bound matrix scale must be positive for "+binding);
		require(!automaticI0 || boundMatrixScale>0,
			"Automatic raw i0 must be bound before validation for "+binding);
		require(!automaticI0 || boundMatrixScale==matrixScale,
			"Automatic raw i0 requires the same matrix scale in every binding: "+binding);
		final double effectiveI0;
		if(requiresView()){
			if(view==null){throw new IllegalArgumentException("Raw per-view formula guard requires a LONG or SHORT binding for "+binding);}
			effectiveI0=(view==View.LONG ? longI0 : shortI0);
		}else{effectiveI0=i0;}
		final boolean matches=effectiveI0==matrixScale;
		if(!matches && !explicitOverride){
			throw new IllegalArgumentException("Raw alignment-score formula i0 must equal bound matrix scale for "+binding+
				": matrix scale="+matrixScale+" formula i0="+effectiveI0+
				"; set tailrawscorei0override=t only for an explicitly reviewed mismatch");
		}
	}

	public static void main(String[] args){
		try{
			final Builder builder=new Builder();
			for(String arg:args){
				final int split=arg.indexOf('=');
				require(split>0 && split<arg.length()-1,"Formula arguments must be flag=value");
				require(builder.parse(arg.substring(0,split),arg.substring(split+1)),"Unknown formula option: "+arg);
			}
			final FamilyBonusFormula formula=builder.build();
			require(formula!=null,"Formula validation requires a form");
			System.out.print(formula.recipe());
			System.out.println("FORMULA_PARAMETER_PASS");
		}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}

	public float score(FamilyTailCascade.Assignment assignment){
		return score(assignment,null);
	}
	public float score(FamilyTailCascade.Assignment assignment,FamilyMatchInformation.View view){
		return score(assignment,view,-1);
	}
	/** Score one accepted assignment for a candidate ORF. candidateAaLength is the
	 * translated CDS length in amino acids, excluding a real terminal stop codon. */
	public float score(FamilyTailCascade.Assignment assignment,FamilyMatchInformation.View view,int candidateAaLength){
		require(assignment!=null,"Formula scoring requires an assignment");
		if(!assignment.assigned()){return 0;}
		final double value;
		switch(form){
			case AXPLUSB:
				value=a*ratio(assignment)+b;
				break;
			case BLOSUM_RATIO:
				final double g=lengthMode.weight(view,candidateAaLength);
				value=mix(outer,shortOuter,g)*(Math.min(cap,mix(scale,shortScale,g)*ratio(assignment))+
					mix(a,shortA,g)+mix(b,shortB,g)*overlap(assignment));
				break;
			case SQRT_RATIO:
				value=Math.min(cap,scale*Math.sqrt(ratio(assignment,view)))+a+b*overlap(assignment);
				break;
			case X_Y_MIN_Z_RATIO:
				value=x+y*Math.min(z,ratio(assignment));
				break;
			default:
				throw new AssertionError("Unhandled bonus formula: "+form);
		}
		require(Double.isFinite(value) && value<=Float.MAX_VALUE && (form.signed || value>=0) &&
			(!form.signed || value>=-Float.MAX_VALUE),
			form.signed ? "Signed formula bonus must be finite and within the caller float range" :
			"Formula bonus must be finite, nonnegative and within the caller float range");
		return (float)value;
	}

	private double ratio(FamilyTailCascade.Assignment assignment){
		final double value=signal(assignment)/i0;
		require(Double.isFinite(value) && value>=0,"Formula signal divided by i0 must be finite and nonnegative");
		return value;
	}
	private double ratio(FamilyTailCascade.Assignment assignment,FamilyMatchInformation.View view){
		final double value=signal(assignment)/viewI0(view);
		require(Double.isFinite(value) && value>=0,"Formula signal divided by view i0 must be finite and nonnegative");
		return value;
	}

	private double signal(FamilyTailCascade.Assignment assignment){
		switch(signal){
			case BLOSUM_RAW:
				require(assignment.rawScore>0,"Raw-BLOSUM scoring requires a positive accepted alignment score");
				return assignment.rawScore;
			case EXACT_INFORMATION:
				require(assignment.informationMeasured,"Exact-information scoring requires measured information evidence");
				return assignment.informationRatio;
			case ALIGNMENT_BLOSUM:
				require(assignment.rawScore>0 && assignment.denominatorSelf>0,
					"Alignment-BLOSUM scoring requires a positive accepted alignment score and self score");
				return assignment.rawScore/(double)assignment.denominatorSelf;
			default:
				throw new AssertionError("Unhandled bonus signal: "+signal);
		}
	}
	private double viewI0(FamilyMatchInformation.View view){
		if(!requiresView()){return i0;}
		require(view!=null,"Per-view sqrt scoring requires an explicit LONG or SHORT family-set view");
		return view==FamilyMatchInformation.View.LONG ? longI0 : shortI0;
	}
	private double overlap(FamilyTailCascade.Assignment assignment){
		if(b==0 && shortB==0){return 0;}
		require(assignment.querySpan>0,"Overlap bonus requires a measured aligned query span");
		return assignment.querySpan;
	}
	private static double mix(double longValue,double shortValue,double shortWeight){
		return longValue+(shortValue-longValue)*shortWeight;
	}

	public ByteBuilder recipe(){
		final ByteBuilder out=new ByteBuilder("form\t").append(form.label).nl();
		out.append("signal\t").append(signal.label).nl();
		if(form==Form.AXPLUSB){
			out.append("i0\t").append(Double.toString(i0)).nl();
			out.append("a\t").append(Double.toString(a)).nl().append("b\t").append(Double.toString(b)).nl();
		}else if(form==Form.BLOSUM_RATIO){
			out.append("i0\t").append(Double.toString(i0)).nl();
			out.append("scale\t").append(Double.toString(scale)).nl().append("cap\t").append(Double.toString(cap)).nl()
				.append("a\t").append(Double.toString(a)).nl().append("b\t").append(Double.toString(b)).nl();
			if(outer!=1){out.append("outer\t").append(Double.toString(outer)).nl();}
			if(lengthMode!=LengthMode.OFF){
				out.append("formulalengthmode\t").append(lengthMode.label).nl()
					.append("shortscale\t").append(Double.toString(shortScale)).nl()
					.append("shorta\t").append(Double.toString(shortA)).nl()
					.append("shortb\t").append(Double.toString(shortB)).nl()
					.append("shortouter\t").append(Double.toString(shortOuter)).nl();
			}
		}else if(form==Form.SQRT_RATIO){
			if(requiresView()){
				out.append("longi0\t").append(Double.toString(longI0)).nl()
					.append("shorti0\t").append(Double.toString(shortI0)).nl();
			}else{out.append("i0\t").append(Double.toString(i0)).nl();}
			out.append("scale\t").append(Double.toString(scale)).nl().append("cap\t").append(Double.toString(cap)).nl()
				.append("a\t").append(Double.toString(a)).nl().append("b\t").append(Double.toString(b)).nl();
		}else if(form==Form.X_Y_MIN_Z_RATIO){
			out.append("i0\t").append(Double.toString(i0)).nl();
			out.append("x\t").append(Double.toString(x)).nl().append("y\t").append(Double.toString(y)).nl()
				.append("z\t").append(Double.toString(z)).nl();
		}else{throw new AssertionError("Unhandled bonus formula: "+form);}
		return out;
	}

	public static final class Builder{
		public boolean parse(String key,String value){
			if(key==null){throw new IllegalArgumentException("Formula option key is required");}
			key=key.toLowerCase(Locale.ROOT);
			if(key.startsWith("tail")){key=key.substring(4);}
			if(!isKey(key)){return false;}
			require(value!=null && !value.isEmpty(),"Formula option "+key+" requires an explicit value");
			require(seen.add(key),"Duplicate formula option: "+key);
			if(key.equals("form")){form=Form.parse(value);}
			else if(key.equals("signal")){signal=Signal.parse(value);}
			else if(key.equals("a")){a=finite(value,key);}
			else if(key.equals("b")){b=finite(value,key);}
			else if(key.equals("x")){x=nonnegative(value,key);}
			else if(key.equals("y")){y=nonnegative(value,key);}
			else if(key.equals("z")){z=positive(value,key);}
			else if(key.equals("scale")){scale=nonnegative(value,key);}
			else if(key.equals("cap")){cap=positive(value,key);}
			else if(key.equals("outer")){outer=positive(value,key);}
			else if(key.equals("shorta")){shortA=finite(value,key);}
			else if(key.equals("shortb")){shortB=finite(value,key);}
			else if(key.equals("shortscale")){shortScale=nonnegative(value,key);}
			else if(key.equals("shortouter")){shortOuter=positive(value,key);}
			else if(key.equals("formulalengthmode")){lengthMode=LengthMode.parse(value);}
			else if(key.equals("i0")){i0=positive(value,key);}
			else if(key.equals("longi0")){longI0=positive(value,key);}
			else if(key.equals("shorti0")){shortI0=positive(value,key);}
			else{throw new AssertionError("Unhandled formula option: "+key);}
			return true;
		}
		public FamilyBonusFormula build(){
			if(seen.isEmpty()){return null;}
			require(form!=null,"Formula constants require a form");
			if(signal==null){signal=Signal.BLOSUM_RAW;}
			if(!seen.contains("i0")){i0=1;}
			if(form==Form.AXPLUSB){
				only("form","signal","a","b","i0");
				require(seen.contains("a") && seen.contains("b"),"AXplusB requires a and b");
				require(a>=0 && b>=0,"AXplusB requires nonnegative a and b");
			}else if(form==Form.BLOSUM_RATIO){
				only("form","signal","scale","cap","outer","shortscale","shortouter",
					"i0","a","b","shorta","shortb","formulalengthmode");
				require(seen.contains("scale") && seen.contains("cap"),"BLOSUMratio requires scale and cap");
				if(lengthMode==null){lengthMode=LengthMode.OFF;}
				final boolean shortEndpoint=seen.contains("shortscale") || seen.contains("shortouter") ||
					seen.contains("shorta") || seen.contains("shortb");
				if(lengthMode==LengthMode.OFF){
					require(!shortEndpoint,"Short formula endpoint requires tailformulalengthmode");
				}else{
					require(seen.contains("shortscale") && seen.contains("shortouter") &&
							seen.contains("shorta") && seen.contains("shortb"),
						"Length-mode BLOSUMratio requires shortscale, shortouter, shorta and shortb");
				}
			}else if(form==Form.SQRT_RATIO){
				only("form","signal","scale","cap","i0","longi0","shorti0","a","b");
				final boolean longSet=seen.contains("longi0"), shortSet=seen.contains("shorti0");
				require(!(seen.contains("i0") && (longSet || shortSet)) &&
					(signal==Signal.BLOSUM_RAW || seen.contains("i0") || longSet || shortSet),
					"sqrtRatio requires either i0 or longi0+shorti0");
				require(longSet==shortSet,"sqrtRatio per-view scoring requires both longi0 and shorti0");
				if(!seen.contains("scale")){scale=400;}
				if(!seen.contains("cap")){cap=400;}
			}else if(form==Form.X_Y_MIN_Z_RATIO){
				only("form","signal","x","y","z","i0");
				require(seen.contains("x") && seen.contains("y") && seen.contains("z"),
					"X+Y*min(Z,ratio) requires x, y and z");
			}else{throw new AssertionError("Unhandled bonus formula: "+form);}
			final boolean automatic=signal==Signal.BLOSUM_RAW && !seen.contains("i0") && !seen.contains("longi0");
			if(lengthMode==null){lengthMode=LengthMode.OFF;}
			return new FamilyBonusFormula(form,signal,lengthMode,a,b,x,y,z,scale,cap,outer,
				shortA,shortB,shortScale,shortOuter,i0,longI0,shortI0,automatic,0);
		}
		private void only(String... allowed){
			for(String key:seen){
				boolean ok=false;for(String a:allowed){ok|=key.equals(a);}
				require(ok,"Formula option "+key+" does not apply to "+form.label);
			}
		}
		private final HashSet<String> seen=new HashSet<String>();
		private Form form;
		private Signal signal;
		private double a,b,x,y,z,scale,cap,outer=1,shortA,shortB,shortScale,shortOuter,i0,longI0,shortI0;
		private LengthMode lengthMode;
	}

	private enum Form{
		AXPLUSB("AXplusB",false), BLOSUM_RATIO("BLOSUMratio",true), SQRT_RATIO("sqrtRatio",true),
		X_Y_MIN_Z_RATIO("X+Y*min(Z,ratio)",false);
		Form(String label_,boolean signed_){label=label_;signed=signed_;}
		static Form parse(String value){
			if("AXplusB".equalsIgnoreCase(value)){return AXPLUSB;}
			if("BLOSUMratio".equalsIgnoreCase(value)){return BLOSUM_RATIO;}
			if("sqrtRatio".equalsIgnoreCase(value)){return SQRT_RATIO;}
			if("cappedRatio".equalsIgnoreCase(value) || "X+Y*min(Z,ratio)".equalsIgnoreCase(value)){
				return X_Y_MIN_Z_RATIO;
			}
			throw new IllegalArgumentException("Unknown bonus formula: "+value);
		}
		final String label;
		final boolean signed;
	}

	private enum LengthMode{
		OFF("off"), VIEW_STEP("familyset_view_step"),
		ORF_STEP("candidate_orf_step"), ORF_RAMP("candidate_orf_ramp");
		LengthMode(String label_){label=label_;}
		/** Short endpoint is length<=99 amino acids; long endpoint is length>=100
		 * for ORF_STEP, and length>=150 for ORF_RAMP. */
		double weight(View view,int candidateAaLength){
			switch(this){
				case OFF:
					return 0;
				case VIEW_STEP:
					require(view!=null,"Formula view-step length mode requires a SHORT or LONG binding");
					return view==View.SHORT ? 1 : 0;
				case ORF_STEP:
					require(candidateAaLength>0,"Formula ORF-step length mode requires a positive candidate amino-acid length");
					return candidateAaLength<=99 ? 1 : 0;
				case ORF_RAMP:
					require(candidateAaLength>0,"Formula ORF-ramp length mode requires a positive candidate amino-acid length");
					return candidateAaLength<=99 ? 1 : candidateAaLength>=150 ? 0 : (150-candidateAaLength)/51.0;
				default:
					throw new AssertionError("Unhandled formula length mode: "+this);
			}
		}
		static LengthMode parse(String value){
			if(value==null){throw new IllegalArgumentException("Formula length mode requires a value");}
			if(value.equalsIgnoreCase("off") || value.equalsIgnoreCase("none")){return OFF;}
			if(value.equalsIgnoreCase("viewstep") || value.equalsIgnoreCase("view_step") ||
					value.equalsIgnoreCase("familyset_view_step")){
				return VIEW_STEP;
			}
			if(value.equalsIgnoreCase("orfstep") || value.equalsIgnoreCase("orf_step") ||
					value.equalsIgnoreCase("candidate_orf_step")){
				return ORF_STEP;
			}
			if(value.equalsIgnoreCase("orframp") || value.equalsIgnoreCase("orf_ramp") ||
					value.equalsIgnoreCase("candidate_orf_ramp")){
				return ORF_RAMP;
			}
			throw new IllegalArgumentException("Unknown formula length mode: "+value);
		}
		final String label;
	}

	private enum Signal{
		BLOSUM_RAW("blosum_raw",false),
		EXACT_INFORMATION("exact_information",true),
		ALIGNMENT_BLOSUM("alignment_blosum",false);
		Signal(String label_,boolean information_){label=label_;information=information_;}
		static Signal parse(String value){
			if(BLOSUM_RAW.label.equalsIgnoreCase(value)){return BLOSUM_RAW;}
			if(EXACT_INFORMATION.label.equalsIgnoreCase(value)){return EXACT_INFORMATION;}
			if(ALIGNMENT_BLOSUM.label.equalsIgnoreCase(value)){return ALIGNMENT_BLOSUM;}
			throw new IllegalArgumentException("Unknown bonus signal: "+value);
		}
		final String label;
		final boolean information;
	}

	private static boolean isKey(String key){
		return key.equals("form") || key.equals("signal") || key.equals("a") || key.equals("b") ||
			key.equals("x") || key.equals("y") || key.equals("z") || key.equals("scale") ||
			key.equals("cap") || key.equals("outer") || key.equals("i0") ||
			key.equals("longi0") || key.equals("shorti0") || key.equals("shorta") ||
			key.equals("shortb") || key.equals("shortscale") || key.equals("shortouter") ||
			key.equals("formulalengthmode");
	}
	private static double finite(String text,String name){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value),name+" must be finite");
		return value;
	}
	private static double nonnegative(String text,String name){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>=0,name+" must be finite and nonnegative");
		return value;
	}
	private static double positive(String text,String name){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>0,name+" must be finite and positive");
		return value;
	}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}

	private final Form form;
	private final Signal signal;
	private final LengthMode lengthMode;
	private final double a,b,x,y,z,scale,cap,outer,i0,longI0,shortI0;
	private final double shortA,shortB,shortScale,shortOuter;
	private final boolean automaticI0;
	private final int boundMatrixScale;
}
