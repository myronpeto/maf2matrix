package org.panda.resource.tcga;

/**
 * This class encodes a single mutation.
 *
 * @author Ozgun Babur
 */
public class MutTuple
{
	public String type;
	public String value;

	public MutTuple(String type, String value)
	{
		this.type = type;
		this.value = value;
	}

	public boolean isDeleterious()
	{
		return value.contains("*") || value.contains("fs") || type.equals("Nonsense");
	}
}
