package test.map;

import java.lang.reflect.Method;

/** Independent signed-value merge cases shared by the ten numeric map variants.
 * @author Shinobu
 */
public final class MapSetToMaxTest{

	/** Reflection keeps the correctness cases identical without changing production interfaces. */
	public static void main(String[] args) throws ReflectiveOperationException{
		check(args.length==1, "Expected expect=baseline|fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation");
		for(String name : new String[]{"IntHashMap", "IntHashMap2", "IntHashMap3",
			"IntLongHashMap2", "IntLongHashMap3", "LongIntMap", "LongIntMapX",
			"LongLongHashMap2", "ObjectIntMap", "ObjectDoubleMap"}){
			checkMap(name, fixed);
		}
	}

	/** Checks missing, overlapping, destination-only and self-merge cases with literal expectations. */
	private static void checkMap(final String name, final boolean fixed) throws ReflectiveOperationException{
		final Class<?> type=Class.forName("map."+name);
		final boolean objects=name.startsWith("Object");
		final Class<?> keyType=objects ? Object.class : name.startsWith("Int") ? int.class : long.class;
		final Class<?> valueType=name.contains("Double") ? double.class :
			(name.startsWith("IntLong") || name.startsWith("LongLong")) ? long.class : int.class;
		final Object source=objects ? type.getConstructor(Class.class).newInstance(String.class) :
			type.getConstructor(int.class).newInstance(64);
		final Object dest=objects ? type.getConstructor(Class.class).newInstance(String.class) :
			type.getConstructor(int.class).newInstance(64);
		final Method put=type.getMethod("put", keyType, valueType), get=type.getMethod("get", keyType);
		final Method merge=type.getMethod("setToMax", type), size=type.getMethod("size");
		final long[] incoming={-5, -5, -5, 7, 0};
		for(int i=0; i<incoming.length; i++){put.invoke(source, key(i, keyType), value(incoming[i], valueType));}
		put.invoke(dest, key(1, keyType), value(-9, valueType));
		put.invoke(dest, key(2, keyType), value(3, valueType));
		put.invoke(dest, key(5, keyType), value(-12, valueType));
		merge.invoke(dest, source);
		final long[] expected={fixed ? -5 : -1, -5, 3, 7, 0, -12};
		check(((Number)size.invoke(dest)).intValue()==6, name+" merged size");
		for(int i=0; i<expected.length; i++){
			check(((Number)get.invoke(dest, key(i, keyType))).longValue()==expected[i], name+" destination key "+i);
		}
		for(int i=0; i<incoming.length; i++){
			check(((Number)get.invoke(source, key(i, keyType))).longValue()==incoming[i], name+" source changed at "+i);
		}
		merge.invoke(dest, dest);
		for(int i=0; i<expected.length; i++){
			check(((Number)get.invoke(dest, key(i, keyType))).longValue()==expected[i], name+" self-merge key "+i);
		}
		if(!objects){
			final Object sentinel=type.getMethod("invalid").invoke(dest);
			put.invoke(source, sentinel, value(-5, valueType));
			merge.invoke(dest, source);
			check(((Number)get.invoke(dest, sentinel)).longValue()==(fixed ? -5 : -1), name+" destination sentinel key");
		}
		checkGrowth(type, objects, keyType, valueType, put, get, merge, size, fixed);
		System.out.println("PASS\t"+name+"\tfixed="+fixed);
	}

	/** Forces destination growth while importing negative values; then verifies self-merge. */
	private static void checkGrowth(final Class<?> type, final boolean objects,
		final Class<?> keyType, final Class<?> valueType, final Method put, final Method get,
		final Method merge, final Method size, final boolean fixed) throws ReflectiveOperationException{
		final Object source=objects ? type.getConstructor(int.class, Class.class).newInstance(64, String.class) :
			type.getConstructor(int.class).newInstance(64);
		final Object dest=objects ? type.getConstructor(int.class, Class.class).newInstance(1, String.class) :
			type.getConstructor(int.class).newInstance(1);
		for(int i=0; i<128; i++){put.invoke(source, key(i, keyType), value(-1L-i, valueType));}
		merge.invoke(dest, source);
		check(((Number)size.invoke(dest)).intValue()==128, type.getSimpleName()+" growing destination size");
		merge.invoke(dest, dest);
		for(int i=0; i<128; i++){
			final long expected=fixed ? -1L-i : -1;
			check(((Number)get.invoke(dest, key(i, keyType))).longValue()==expected,
				type.getSimpleName()+" growth/self-merge key "+i);
			check(((Number)get.invoke(source, key(i, keyType))).longValue()==-1L-i,
				type.getSimpleName()+" growing merge changed source key "+i);
		}
	}

	private static Object key(final int index, final Class<?> type){
		if(type==Object.class){return "key"+index;}
		if(type==int.class){return Integer.valueOf(17+index);}
		return Long.valueOf(17+index);
	}

	private static Object value(final long value, final Class<?> type){
		if(type==double.class){return Double.valueOf(value);}
		if(type==long.class){return Long.valueOf(value);}
		return Integer.valueOf((int)value);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
