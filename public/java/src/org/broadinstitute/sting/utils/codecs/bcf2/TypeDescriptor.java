/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.bcf2;

/**
 * Simple BCF decoder
 * @author Mark DePristo
 * @since 5/3/12
 */
public class TypeDescriptor {
    public static final int OVERFLOW_ELEMENT_MARKER = 15;
    public static final int MAX_INLINE_ELEMENTS = 14;

    public final static BCFType[] INTEGER_TYPES_BY_SIZE = new BCFType[3];
    public final static BCFType[] DICTIONARY_TYPES_BY_SIZE = INTEGER_TYPES_BY_SIZE;
    private final static BCFType[] lookup = BCFType.values();

    static {
        INTEGER_TYPES_BY_SIZE[0] = BCFType.INT8;
        INTEGER_TYPES_BY_SIZE[1] = BCFType.INT16;
        INTEGER_TYPES_BY_SIZE[2] = BCFType.INT32;
    }

    public final static byte encodeTypeDescriptor(final int nElements, final BCFType type ) {
        int encodeSize = Math.min(nElements, OVERFLOW_ELEMENT_MARKER);
        byte typeByte = (byte)((0x0F & encodeSize) << 4 | (type.getID() & 0x0F));
        return typeByte;
    }

    public final static int decodeSize(final byte typeDescriptor) {
        return (0xF0 & typeDescriptor) >> 4;
    }

    public final static int decodeTypeID(final byte typeDescriptor) {
        return typeDescriptor & 0x0F;
    }

    public final static BCFType decodeType(final byte typeDescriptor) {
        return lookup[decodeTypeID(typeDescriptor)];
    }

    public final static boolean sizeIsOverflow(final byte typeDescriptor) {
        return decodeSize(typeDescriptor) == OVERFLOW_ELEMENT_MARKER;
    }

    public final static boolean willOverflow(final long nElements) {
        return nElements > MAX_INLINE_ELEMENTS;
    }
}
