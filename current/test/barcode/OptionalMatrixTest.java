package test.barcode;

import java.lang.reflect.Field;
import java.util.Locale;

import barcode.PCRMatrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks optional matrix selection in a fresh JVM without linking private classes.
 * The load case inspects the factory's selected types and exercises its parsing
 * hooks and public Hamming-distance constructor. Optional constructors run in
 * separate cases so an existing license rejection cannot masquerade as a load check.
 * @author Brian Bushnell, Shinobu
 */
public final class OptionalMatrixTest{

	/** Requires expected=absent|present and case=load|prob|tile. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, OptionalMatrixTest.class, false);
		try{
			String expected=null, name=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT);
				final String value=pair.length>1 ? pair[1] : null;
				if(key.equals("expected")){expected=value;}
				else if(key.equals("case")){name=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			if(!"absent".equals(expected) && !"present".equals(expected)){
				throw new IllegalArgumentException("Set expected=absent or expected=present.");
			}
			if(!"load".equals(name) && !"prob".equals(name) && !"tile".equals(name)){
				throw new IllegalArgumentException("Set case=load, case=prob or case=tile.");
			}
			final boolean present=expected.equals("present");
			if(name.equals("load")){checkLoading(present);}
			else{
				final int type=name.equals("prob") ? PCRMatrix.PROB_TYPE : PCRMatrix.TILE_TYPE;
				final PCRMatrix matrix=PCRMatrix.create(type, 4, 4, '+', false);
				if(!present){throw new AssertionError("An absent implementation must not construct a usable matrix.");}
				final String className="barcode.prob."+(name.equals("prob") ? "PCRMatrixProb" : "PCRMatrixTile");
				checkName(matrix.getClass(), className);
				System.out.println("OPTIONAL_MATRIX case="+name+" class="+matrix.getClass().getName());
			}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each process has one explicit classpath and configuration. */
	private OptionalMatrixTest(){}

	/** Verifies all three selected classes rather than merely finding optional files. */
	private static void checkLoading(final boolean present) throws Exception{
		final String prefix=present ? "barcode.prob." : "barcode.stub.";
		final String[] fields={"pcrMatrixProbAbstractClass", "pcrMatrixProbClass", "pcrMatrixTileClass"};
		final String[] names={"PCRMatrixProbAbstract", "PCRMatrixProb", "PCRMatrixTile"};
		for(int i=0; i<fields.length; i++){
			final Field field=PCRMatrix.class.getDeclaredField(fields[i]);
			field.setAccessible(true);
			checkName((Class<?>)field.get(null), prefix+names[i]);
		}
		final boolean parsed=PCRMatrix.parseStatic("accuracy=0.9", "accuracy", "0.9");
		if(parsed!=present){throw new AssertionError("Optional parser selected incorrectly: expected="+present+", parsed="+parsed);}
		PCRMatrix.postParseStatic();
		final PCRMatrix matrix=PCRMatrix.create(PCRMatrix.HDIST_TYPE, 4, 4, '+', false);
		checkName(matrix.getClass(), "barcode.PCRMatrixHDist");
		System.out.println("OPTIONAL_MATRIX case=load expected="+(present ? "present" : "absent")+" selected="+prefix+" checks=5");
	}

	/** Exact class identity is the factory-selection oracle; a subclass is not equivalent here. */
	private static void checkName(final Class<?> actual, final String expected){
		if(!actual.getName().equals(expected)){
			throw new AssertionError("Expected class "+expected+", selected "+actual.getName());
		}
	}
}
