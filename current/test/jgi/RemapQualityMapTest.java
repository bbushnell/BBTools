package test.jgi;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Locale;

import jgi.RemapQuality;
import parse.Parse;
import stream.SamLine;

/**
 * One-case RemapQuality constructor map fixture for signed-byte checks.
 * @author Jean
 */
public final class RemapQualityMapTest{

	/** Requires in=, out=, phase=baseline|fixed, and case=. */
	public static void main(final String[] args) throws Exception{
		String in=null, out=null, phase=null, caseName=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("in")){in=pair[1];}
			else if(key.equals("out")){out=pair[1];}
			else if(key.equals("phase")){phase=pair[1];}
			else if(key.equals("case")){caseName=pair[1];}
			else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase),
			"Set phase=baseline|fixed.");
		require(caseName!=null && !caseName.isEmpty(), "Set case=<name>.");
		require(in!=null && !in.isEmpty(), "Set in=<existing tiny FASTQ>.");
		require(out!=null && !out.isEmpty(), "Set out=<case TSV>.");

		final Path input=Paths.get(in).toAbsolutePath().normalize();
		final Path output=Paths.get(out).toAbsolutePath().normalize();
		require(Files.isRegularFile(input), "Input FASTQ does not exist: "+input);
		new RemapQualityMapTest(input, output, phase, caseName).run();
	}

	/** Records the selected literal case. */
	private void run() throws Exception{
		final Case expected=findCase(caseName);
		final boolean setFromOk=SamLine.SET_FROM_OK;
		Throwable failure=null;
		RemapQuality tool=null;
		boolean restored=false;

		try{
			try{
				tool=new RemapQuality(constructorArgs(input, expected.mapString));
			}catch(final Throwable caught){
				failure=caught;
			}
		}finally{
			SamLine.SET_FROM_OK=setFromOk;
			restored=(SamLine.SET_FROM_OK==setFromOk);
		}

		final Class<? extends Throwable> expectedException=
			"fixed".equals(phase) ? expected.fixedException : null;
		final boolean expectedRejection=(expectedException!=null);
		final boolean actualRejection=(failure!=null);
		final int[] expectedValues=expectedRejection ? null : expected.baselineValues;
		final int[] actualValues=(tool==null ? null : values(tool, expected.indexes));
		final boolean exceptionMatches=expectedRejection ?
			failure!=null && failure.getClass()==expectedException : failure==null;
		final boolean valuesMatch=expectedRejection ||
			Arrays.equals(expectedValues, actualValues);
		final boolean match=exceptionMatches && valuesMatch && restored;

		final String text=HEADER+'\n'+
			expected.name+'\t'+phase+'\t'+operation(expected.mapString)+'\t'+
			clean(expected.mapString==null ? "" : expected.mapString)+'\t'+input+'\t'+
			expectedRejection+'\t'+actualRejection+'\t'+
			Arrays.toString(expected.indexes)+'\t'+
			(expectedValues==null ? "unavailable" : Arrays.toString(expectedValues))+'\t'+
			(actualValues==null ? "unavailable" : Arrays.toString(actualValues))+'\t'+
			exceptionName(expectedException)+'\t'+exceptionName(failure)+'\t'+
			clean(failure==null ? "" : String.valueOf(failure.getMessage()))+'\t'+
			restored+'\t'+match+'\n';
		Files.write(output, text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		if(!match){
			throw new AssertionError("Unexpected RemapQuality constructor observation; "+
				"inspect "+output);
		}
	}

	/** Selects one of the twelve literal map cases. */
	private static Case findCase(final String name){
		for(final Case row : CASES){
			if(row.name.equals(name)){return row;}
		}
		throw new IllegalArgumentException("Unknown case: "+name);
	}

	/** Builds the common constructor arguments plus one optional literal map. */
	private static String[] constructorArgs(final Path input, final String mapString){
		final String[] common={"in="+input, "out=null", "threads=1", "simd=f"};
		if(mapString==null){return common;}
		return new String[]{common[0], common[1], common[2], common[3],
			"map="+mapString};
	}

	/** Reads the signed byte value at every requested unsigned index. */
	private static int[] values(final RemapQuality tool, final int[] indexes){
		final int[] values=new int[indexes.length];
		for(int i=0; i<indexes.length; i++){
			values[i]=tool.map[indexes[i]];
		}
		return values;
	}

	/** Names the constructor path under test. */
	private static String operation(final String mapString){
		return mapString==null ? "constructor-default-map" : "constructor-custom-map";
	}

	/** Makes one value fit in a TSV cell. */
	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}

	/** Records absent exceptions explicitly. */
	private static String exceptionName(final Class<? extends Throwable> clazz){
		return clazz==null ? "none" : clazz.getName();
	}

	/** Records absent exceptions explicitly. */
	private static String exceptionName(final Throwable error){
		return error==null ? "none" : error.getClass().getName();
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	private RemapQualityMapTest(final Path input_, final Path output_,
			final String phase_, final String caseName_){
		input=input_;
		output=output_;
		phase=phase_;
		caseName=caseName_;
	}

	private final Path input;
	private final Path output;
	private final String phase;
	private final String caseName;

	private static final Class<? extends Throwable> NFE=
		NumberFormatException.class;

	private static final Case[] CASES={
		new Case("default", null,
			new int[]{0, 1, 2, 41, 42, 127, 128, 255},
			new int[]{0, 1, 41, 2, 42, 127, -128, -1}, null),
		new Case("custom_identity", "2,2",
			new int[]{2, 3, 41}, new int[]{2, 3, 41}, null),
		new Case("duplicate_trailing", "2,3;2,4;",
			new int[]{2}, new int[]{4}, null),
		new Case("min_byte", "2,-128",
			new int[]{2}, new int[]{-128}, null),
		new Case("max_byte", "2,+127",
			new int[]{2}, new int[]{127}, null),
		new Case("source_255", "255,7",
			new int[]{255}, new int[]{7}, null),
		new Case("overflow_high", "2,128",
			new int[]{2}, new int[]{-128}, NFE),
		new Case("overflow_low", "2,-129",
			new int[]{2}, new int[]{127}, NFE),
		new Case("wrap_positive", "2,256",
			new int[]{2}, new int[]{0}, NFE),
		new Case("wrap_negative", "2,-256",
			new int[]{2}, new int[]{0}, NFE),
		new Case("int_max", "2,2147483647",
			new int[]{2}, new int[]{-1}, NFE),
		new Case("int_min", "2,-2147483648",
			new int[]{2}, new int[]{0}, NFE)
	};

	private static final String HEADER="case	phase	operation	map	input"+
		"	expectedRejection	actualRejection	indexes"+
		"	expectedValues	actualValues	expectedException"+
		"	actualException	exceptionMessage	setFromOkRestored"+
		"	expectedObservationMatched";

	private static final class Case{
		Case(final String name_, final String mapString_, final int[] indexes_,
				final int[] baselineValues_,
				final Class<? extends Throwable> fixedException_){
			name=name_;
			mapString=mapString_;
			indexes=indexes_;
			baselineValues=baselineValues_;
			fixedException=fixedException_;
		}

		final String name;
		final String mapString;
		final int[] indexes;
		final int[] baselineValues;
		final Class<? extends Throwable> fixedException;
	}
}
