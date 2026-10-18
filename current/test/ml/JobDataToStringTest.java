package test.ml;

import ml.JobData;
import ml.JobResults;
import ml.Sample;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

/**
 * Selectable direct-call fixture for JobData.toString() backing-source diagnostics.
 *
 * @author Sucrose
 */
public final class JobDataToStringTest{

	/** Runs one named JobData.toString() case. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Run one case name");}
		run(args[0]);
		System.err.println("JobDataToStringTest PASS "+args[0]);
	}

	/** Dispatches the requested direct-call case. */
	private static void run(final String name){
		if(name.equals("list")){
			check("list", listJob(), "jD: epoch=3, max=1, len=2, back=true, sort=false");
		}else if(name.equals("set")){
			check("set", setJob(), "jD: epoch=4, max=2, len=3, back=false, sort=false");
		}else if(name.equals("poison")){
			check("poison", JobData.POISON, "jD: epoch=-1, max=-1, len=0, back=false, sort=false");
		}else{
			throw new IllegalArgumentException("Unknown fixture mode: "+name);
		}
	}

	/** Builds a list-backed job whose available list is larger than its processing limit. */
	private static JobData listJob(){
		final ArrayList<Sample> list=new ArrayList<Sample>();
		list.add(sample(0));
		list.add(sample(1));
		return new JobData(null, queue(), 3, 1, 0.5,
			true, 1.0f, false, false, list, null, null, -1, 1);
	}

	/** Builds an array-backed job whose backing set is larger than its processing limit. */
	private static JobData setJob(){
		final Sample[] set=new Sample[]{sample(0), sample(1), sample(2)};
		return new JobData(null, queue(), 4, 2, 0.5,
			false, 1.0f, false, false, null, set, null, -1, 1);
	}

	/** Creates a real queue so normal jobs are direct constructor calls. */
	private static ArrayBlockingQueue<JobResults> queue(){
		return new ArrayBlockingQueue<JobResults>(1);
	}

	/** Builds a minimal real Sample; JobData.toString() uses only collection sizes. */
	private static Sample sample(final int id){
		return new Sample(new float[]{id}, new float[]{id&1}, 1, id);
	}

	/** Checks the exact diagnostic string for one backing-source case. */
	private static void check(final String name, final JobData job, final String expected){
		final String actual=job.toString();
		if(!expected.equals(actual)){
			throw new AssertionError(name+" expected '"+expected+"' but observed '"+actual+"'");
		}
	}

	private JobDataToStringTest(){}

}
