package synth;

import java.lang.reflect.Field;

/** Observes constructed pair spans without entering the shredding/copying path.
 * @author Brian Bushnell, Jean
 */
public final class ShredPacBioTest{

	/** Emits five observations for comparison with literal external span oracles.
	 * @param args One in= argument naming a small existing input for format validation
	 * @throws ReflectiveOperationException If the production span field cannot be inspected
	 */
	public static void main(String[] args) throws ReflectiveOperationException{
		require(args.length==1 && args[0].startsWith("in=") && args[0].length()>3,
				"Specify in=<fixture> for ShredPacBio constructor validation");
		require(ShredPacBio.class.desiredAssertionStatus(), "Production assertions must be enabled");
		final Field span=ShredPacBio.class.getDeclaredField("fragLength");
		span.setAccessible(true);
		System.out.println("read_length\tfragment_span");
		observe(span, args[0], 1);
		observe(span, args[0], 500);
		observe(span, args[0], 1073741823);
		observe(span, args[0], 1073741824);
		observe(span, args[0], Integer.MAX_VALUE);
	}

	/** Constructs only; no process method, stream startup or sequence copying occurs. */
	private static void observe(Field span, String input, int readLength) throws IllegalAccessException{
		require(readLength>0, "Constructor requires positive readLength for the stored pair span");
		final ShredPacBio shredder=new ShredPacBio(new String[]{input, "out=null",
				"readlen="+readLength, "minlen=1", "t=1", "simd=f", "bf1=t"});
		final Object value=span.get(shredder);
		require(value instanceof Number, "fragLength must be numeric for the production pair-span comparison");
		System.out.println(readLength+"\t"+((Number)value).longValue());
	}

	/** Keeps fixture preconditions active independently of assertion mode. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Static-only fixture entry point. */
	private ShredPacBioTest(){}
}
