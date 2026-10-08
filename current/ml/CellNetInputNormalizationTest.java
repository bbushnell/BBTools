package ml;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import structures.FloatList;
import structures.ByteBuilder;

/** Exercises explicit input normalization through parsing, copies and serialization.
 * @author Yoimiya
 */
public final class CellNetInputNormalizationTest {

	/** Runs deterministic analytic fixtures; no training or benchmark is performed. */
	public static void main(String[] args){
		assert(args.length==0) : "This analytic fixture accepts no project parameters";
		for(boolean dense : new boolean[]{true,false}){
			final String text=fixture(dense);
			final CellNet net=parse(text);
			check(net);
			final CellNet copied=net.copy(false),denseCopy=net.copyDenseForInference();
			check(copied);
			check(denseCopy);
			final CellNet assigned=parse(text);
			assigned.setFrom(net,false);
			check(assigned);
			CellNet.codingA48Out=true;
			check(parse(new String(net.toBytes().toBytes(),StandardCharsets.UTF_8)));
			final float[] mean={16777216f,1f},scale={1e9f,0.5f};
			net.setInputNormalization(mean,scale);
			mean[0]=0; scale[1]=9;
			check(net);
			net.setInputNormalization(null,null);
			if(net.hasInputNormalization()){throw new AssertionError("Identity reset retained preprocessing");}
			check(copied); check(denseCopy); check(assigned);
			net.applyInput(new float[]{1,2});
			if(net.values[0][0]!=1 || net.values[0][1]!=2){throw new AssertionError("Historical input path changed");}
			final String meanHeader=header("#inputmean_a48",16777216f,1f);
			final String scaleHeader=header("#inputinversestd_a48",1e9f,0.5f);
			reject(text.replace(meanHeader,""));
			reject(text.replace(scaleHeader,header("#inputinversestd_a48",1e9f,0f)));
			reject(text.replace(meanHeader,header("#inputmean_a48",16777216f)));
			reject(text.replace(meanHeader,header("#inputmean_a48")));
			reject(text.replace(meanHeader,meanHeader+meanHeader));
			reject(text.replace(meanHeader+scaleHeader,""));
		}
		System.out.println("CellNetInputNormalizationTest PASS");
	}

	/** A large mean and scale exercise the cancellation that prevents safe weight folding. */
	private static void check(CellNet net){
		assert(net!=null) : "Normalization checks require a constructed network";
		final float[] row={16777216f,2f};
		net.applyInput(row);
		if(!net.hasInputNormalization() || net.values[0][0]!=0 || net.values[0][1]!=0.5f){
			throw new AssertionError("Subtract-then-scale input arithmetic changed");
		}
		net.feedForward();
		if(net.getOutput(0)!=3f || row[0]!=16777216f || row[1]!=2f){
			throw new AssertionError("Normalized inference changed or mutated its caller's row");
		}
		final FloatList list=new FloatList(2);
		list.add(row[0]); list.add(row[1]);
		net.applyInput(list);
		net.feedForward();
		if(net.getOutput(0)!=3f){throw new AssertionError("FloatList preprocessing differs from float[]");}
	}

	/** Parses exactly the same text format used by normal CellNet file loading. */
	private static CellNet parse(String text){
		assert(text!=null) : "A parser fixture must supply its serialized network";
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(String line : text.split("\n")){lines.add(line.getBytes(StandardCharsets.UTF_8));}
		return CellNetParser.loadFromLines(lines);
	}

	/** Malformed preprocessing must fail before a network can be used for inference. */
	private static void reject(String text){
		assert(text!=null) : "Negative fixture needs the malformed bytes it is testing";
		try{parse(text);}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed input normalization was silently accepted");
	}

	/** Constructs a two-input linear net whose standardized expected output is exactly three. */
	private static String fixture(boolean dense){
		return "##bbnet\n#version 1\n#concise\n"+(dense ? "#dense\n" : "#sparse\n")+
			"#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 2 1\n"+
			header("#inputmean_a48",16777216f,1f)+header("#inputinversestd_a48",1e9f,0.5f)+
			"##normalization explicit-input\n#coding decimal\n"+
			(dense ? "C3 LINEAR 1 3 4\n" : "I3 0 1\nW3 LINEAR 1 3 4\n");
	}

	/** Encodes fixture constants with the native exact-bit header coding. */
	private static String header(String name,float... values){
		assert(name!=null) : "Fixture headers require their critical field name";
		final ByteBuilder out=new ByteBuilder();
		out.append(name);
		for(float value : values){out.space().appendFloatA48(value);}
		return out.nl().toString();
	}
}
