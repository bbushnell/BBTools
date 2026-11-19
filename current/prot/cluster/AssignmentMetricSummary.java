package prot.cluster;

import java.util.Arrays;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Final-pass accepted member metrics against the last assignment reference.
 * Final consensuses are rebuilt after assignment; these values are not rescored.
 */
public final class AssignmentMetricSummary {

	public AssignmentMetricSummary(int inputs){
		if(inputs<0){throw new IllegalArgumentException("Negative input count: "+inputs);}
		familyByInput=new int[inputs];
		Arrays.fill(familyByInput, -1);
		identities=new double[inputs];
		scoreRatios=new double[inputs];
	}

	public void record(int inputIndex,int familyId,double identity,double scoreRatio){
		if(inputIndex<0 || inputIndex>=familyByInput.length){
			throw new IllegalArgumentException("Input index out of range: "+inputIndex);
		}
		if(familyId<0){throw new IllegalArgumentException("Negative family ID: "+familyId);}
		if(familyByInput[inputIndex]>=0){
			throw new IllegalStateException("Duplicate assignment metric for input "+inputIndex);
		}
		if(!Double.isFinite(identity) || identity<0 || identity>1){
			throw new IllegalArgumentException("Identity must be finite in [0,1]: "+identity);
		}
		if(!Double.isFinite(scoreRatio)){
			throw new IllegalArgumentException("Score ratio must be finite: "+scoreRatio);
		}
		familyByInput[inputIndex]=familyId;
		identities[inputIndex]=identity;
		scoreRatios[inputIndex]=scoreRatio;
		acceptedBeforePrune++;
	}

	public int acceptedBeforePrune(){return acceptedBeforePrune;}

	public int inputs(){return familyByInput.length;}

	public void write(String path,boolean overwrite,int[] finalAssignments){
		if(path==null || path.length()<1){throw new IllegalArgumentException("Metric summary output path is required");}
		if(finalAssignments==null || finalAssignments.length!=familyByInput.length){
			throw new IllegalArgumentException("Final assignment count differs from metric input count");
		}
		final double[] retainedIdentity=new double[acceptedBeforePrune];
		final double[] retainedRatio=new double[acceptedBeforePrune];
		int retained=0;
		for(int input=0; input<familyByInput.length; input++){
			final int finalFamily=finalAssignments[input], recordedFamily=familyByInput[input];
			if(finalFamily>=0 && finalFamily!=recordedFamily){
				throw new IllegalStateException("Retained assignment lacks matching recorded metric for input "+input);
			}
			if(finalFamily>=0){
				retainedIdentity[retained]=identities[input];
				retainedRatio[retained]=scoreRatios[input];
				retained++;
			}
		}
		if(retained<0 || retained>acceptedBeforePrune){
			throw new IllegalStateException("Retained metric count exceeded accepted winners");
		}
		Arrays.sort(retainedIdentity,0,retained);
		Arrays.sort(retainedRatio,0,retained);
		final double identityMean=retained>0 ? mean(retainedIdentity,retained) : 0;
		final double ratioMean=retained>0 ? mean(retainedRatio,retained) : 0;
		final double identityMedian=retained>0 ? median(retainedIdentity,retained) : 0;
		final double ratioMedian=retained>0 ? median(retainedRatio,retained) : 0;
		if(retained>0 && !(Double.isFinite(identityMean) && Double.isFinite(ratioMean) &&
				Double.isFinite(identityMedian) && Double.isFinite(ratioMedian))){
			throw new IllegalStateException("Assignment metric summary produced a nonfinite statistic");
		}
		final ByteBuilder bb=new ByteBuilder();
		bb.append("#measurement\tlast_assignment_reference\n");
		bb.append("#consensus\tfinal consensus is rebuilt after these assignments; metrics were not rescored against that rebuilt consensus\n");
		bb.append("accepted_before_prune\tretained_members\tvalue_available\tidentity_mean\tidentity_median\tscore_ratio_mean\tscore_ratio_median\n");
		bb.append(acceptedBeforePrune).tab().append(retained).tab().append(retained>0 ? 1 : 0).tab();
		if(retained>0){
			appendMetric(bb,identityMean).tab(); appendMetric(bb,identityMedian).tab();
			appendMetric(bb,ratioMean).tab(); appendMetric(bb,ratioMedian);
		}else{bb.append("NA\tNA\tNA\tNA");}
		bb.nl();
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false);
		writer.start();
		writer.print(bb);
		if(writer.poisonAndWait()){throw new IllegalStateException("Assignment metric summary output failed: "+path);}
	}

	private static double median(double[] sorted,int count){
		if(count<=0 || count>sorted.length){
			throw new IllegalArgumentException("Invalid median count "+count+" for "+sorted.length+" values");
		}
		final int mid=count/2;
		return (count&1)==1 ? sorted[mid] : sorted[mid-1]*0.5+sorted[mid]*0.5;
	}

	private static double mean(double[] values,int count){
		if(count<=0 || count>values.length){
			throw new IllegalArgumentException("Invalid mean count "+count+" for "+values.length+" values");
		}
		double mean=0;
		for(int i=0; i<count; i++){mean+=values[i]/count;}
		return mean;
	}

	private static ByteBuilder appendMetric(ByteBuilder bb,double value){
		return bb.appendSlow(value);
	}

	private final int[] familyByInput;
	private final double[] identities,scoreRatios;
	private int acceptedBeforePrune;
}
