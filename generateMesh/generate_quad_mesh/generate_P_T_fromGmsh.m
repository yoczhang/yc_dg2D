function [ P, T ] = generate_P_T_fromGmsh( gmsh_filename )
% generate P, T from 'gmsh .msh' formate
% such as: [ P, T ] = generate_P_T_fromGmsh('2D.msh')

[ node_num, node_dim, element_num, element_order ] = gmsh_size_read_1 ( gmsh_filename );

P=zeros(2, node_num);
T=zeros(3, element_num);

[ node_x, element_node ] = gmsh_data_read_1 ( gmsh_filename, node_dim, node_num, element_order, element_num );
% the node_x is the coordinates of one point, but has (x, y, z) z=0, so P
% is the first two rows of node_x.
% element_node is as the T, but in one column is arranged in clockwise, so
% we need to change it in anticlockwise.

P(1, :)=node_x(1, :);
P(2, :)=node_x(2, :);
T(1, :)=element_node(1, :);
T(2, :)=element_node(3, :);
T(3, :)=element_node(2, :);

return
end

function [ node_x, element_node ] = gmsh_data_read_1 ( gmsh_filename, node_dim, node_num, element_order, element_num )

%*****************************************************************************80
%
%% GMSH_DATA_READ reads data from a GMSH file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 October 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string GMSH_FILENAME, the GMSH filename.
%
%    Input, integer NODE_DIM, the spatial dimension.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Output, real NODE_X(NODE_DIM,NODE_NUM), the node coordinates.
%
%    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
%    the nodes that make up each element.
%
  node_x = zeros ( node_dim, node_num );
  element_node = zeros ( element_order, element_num );
%
%  Open the file.
%
  input = fopen ( gmsh_filename, 'rt' );

  if ( input < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_DATA_READ - Error!\n' );
    fprintf ( 1, '  Could not open the input file "%s".\n', gmsh_filename );
    error ( 'GMSH_DATA_READ - Error!' );
    return
  end

  level = 0;

  while ( 1 )
    text = fgetl ( input );

    if ( text == -1 )
      break
    end

    if ( level == 0 )
      if ( s_begin ( text(1:6), '$Nodes' ) )
        level = 1;
        j = 0;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      level = 2;
    elseif ( level == 2 )
        flag_=0;
        if(s_len_trim(text)==9)
            if(s_begin ( text(1:9), '$EndNodes' ))
                flag_=1;
            end
        end
      if ( flag_ )
        break
      else
        j = j + 1;
        [ value, count ] = sscanf ( text, '%d  %f  %f  %f' );
        indx = value(1);
        node_x(1,j) = value(2);
        if ( 2 < count )
          node_x(2,j) = value(3);
          if ( 3 < count )
            node_x(3,j) = value(4);
          end
        end
      end
    end
  end
%
%  Now read element information.
%
  level = 0;

  while ( 1 )

    text = fgetl ( input );

    if ( text == -1 )
      fprintf ( 'ran out\n' );
      break
    end

    if ( level == 0 )
      if ( s_begin ( text(1:9), '$Elements' ) )
        level = 1;
        j = 0;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      level = 2;
    elseif ( level == 2 )
      if ( s_begin ( text(1:12), '$EndElements' ) )
        break
      else
%         j = j + 1;
        [ value, count ] = sscanf ( text, '%d' );
        if(value(2)==2)
            j=j+1;
            for i = 1 : element_order
              element_node(i,j) = value(5+i);
            end
        end
      end
    end

  end

  fclose ( input );

  return
end


function [ node_num, node_dim, element_num, element_order ] = gmsh_size_read_1 ( gmsh_filename )

%*****************************************************************************80
%
%% GMSH_SIZE_READ reads sizes from a GMSH file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 October 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string GMSH_FILENAME, the GMSH filename.
%
%    Output, integer NODE_NUM, the number of nodes.
%
%    Output, integer NODE_DIM, the spatial dimension.
%
%    Output, integer ELEMENT_NUM, the number of elements.
%
%    Output, integer ELEMENT_ORDER, the order of the elements.
%
  r8_big = 1.0E+30;

  node_num = 0;
  node_dim = 0;
  element_num = 0;
  element_order = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;
%
%  Open the file.
%
  input = fopen ( gmsh_filename, 'rt' );

  if ( input < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_SIZE_READ - Error!\n' );
    fprintf ( 1, '  Could not open the input file "%s".\n', gmsh_filename );
    error ( 'GMSH_SIZE_READ - Error!' );
    return
  end

  level = 0;

  while ( 1 )

    text = fgetl ( input );

    if ( text == -1 )
      break
    end

    if ( level == 0 )
        %length_nodes=s_len_trim(text)
      if ( s_begin ( text(1:6), '$Nodes' ) )
        level = 1;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      node_num = value(1);
      level = 2;
    elseif ( level == 2 )
        %length_endnodes=s_len_trim('$EndNodes')
        flag_=0;
        if(s_len_trim(text)==9)
            if(s_begin ( text(1:9), '$EndNodes' ))
                flag_=1;
            end
        end
      if ( flag_ )
        break
      else
        [ value, count ] = sscanf ( text, '%d  %f  %f  %f' );
        indx = value(1);
        x = value(2);
        x_min = min ( x_min, x );
        x_max = max ( x_max, x );
        if ( 2 < count )
          y = value(3);
          y_min = min ( y_min, y );
          y_max = max ( y_max, y );
          if ( 3 < count )
            z = value(4);
            z_min = min ( z_min, z );
            z_max = max ( z_max, z );
          end
        end
      end
    end

  end
%
%  Make a very simple guess as to the dimensionality of the data.
%
  node_dim = 3;
  if ( z_max == z_min )
    node_dim = 2;
    if ( y_max == y_min )
      node_dim = 1;
    end
  end
%
%  Now read element information.
%
  level = 0;

 % text = fgetl ( input );
 %fprintf ( 1, '  text= "%s".\n', text );
  counting = -1;%
  while ( 1 )

    text = fgetl ( input );

    if ( text == -1 )
      fprintf ( 'ran out\n' );
      break
    end

    if ( level == 0 )
        %length_elements=s_len_trim(text)
      if ( s_begin ( text(1:9), '$Elements' ) )
        level = 1;
      end
    elseif ( level == 1 )
      [ value, count ] = sscanf ( text, '%d' );
      element_num = value(1);
      level = 2;
    elseif ( level == 2 )
        %length_endelements=s_len_trim(text)
      if ( s_begin ( text(1:12), '$EndElements' ) )
        break
      else
          counting=counting+1;
        [ value, count ] = sscanf ( text, '%d' );
        if(value(2)==2)
%             value
%             count
            element_order = count - 5;
            element_num=element_num-counting;
            break
        end
%         value
%         count
%         element_order = count - 5;
%         break
      end
    end

  end

  fclose ( input );

  return
end


function value = s_begin ( s1, s2 )

%*****************************************************************************80
%
%% S_BEGIN is TRUE if one string matches the beginning of the other.
%
%  Discussion:
%
%    The strings are compared, ignoring blanks and capitalization.
%
%  Example:
%
%     S1              S2      S_BEGIN
%
%    'Bob'          'BOB'     TRUE
%    '  B  o b '    ' bo b'   TRUE
%    'Bob'          'Bobby'   TRUE
%    'Bobo'         'Bobb'    FALSE
%    ' '            'Bob'     FALSE    (Do not allow a blank to match
%                                       anything but another blank string.)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 January 2006
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character S1(*), S2(*), the strings to be compared.
%
%    Output, logical S_BEGIN, is TRUE if the strings match up to
%    the end of the shorter string, ignoring case.
%
  len1 = s_len_trim ( s1 );
  len2 = s_len_trim ( s2 );
%
%  If either string is blank, then both must be blank to match.
%  Otherwise, a blank string matches anything, which is not
%  what most people want.
%
  if ( len1 == 0 || len2 == 0 )

    if ( len1 == 0 && len2 == 0 )
      value = 1;
    else
      value = 0;
    end

    return

  end

  i1 = 0;
  i2 = 0;
%
%  Find the next nonblank in S1.
%
  while ( 1 )

    while ( 1 )

      i1 = i1 + 1;

      if ( len1 < i1 )
        value = 1;
        return
      end

      if ( s1(i1) ~= ' ' )
        break
      end

    end
%
%  Find the next nonblank in S2.
%
    while ( 1 )

      i2 = i2 + 1;

      if ( len2 < i2 )
        value = 1;
        return
      end

      if ( s2(i2) ~= ' ' )
        break
      end

    end
%
%  If the characters match, get the next pair.
%
    if ( ~ch_eqi ( s1(i1), s2(i2) ) )
      break
    end

  end

  value = 0;

  return
end


function len = s_len_trim ( s )

%*****************************************************************************80
%
%% S_LEN_TRIM returns the length of a character string to the last nonblank.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 June 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string S, the string to be measured.
%
%    Output, integer LEN, the length of the string up to the last nonblank.
%
  len = length ( s );

  while ( 0 < len )
    if ( s(len) ~= ' ' )
      return
    end
    len = len - 1;
  end

  return
end

function c = ch_cap ( c )
%*****************************************************************************80
%
%% CH_CAP capitalizes a single character.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 April 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character C, the character to capitalize.
%
%    Output, character C, the capitalized character.
%
if ( 'a' <= c && c <= 'z' )
    c = c + 'A' - 'a';
end

c = char ( c );

return
end

function truefalse = ch_eqi ( c1, c2 )

%*****************************************************************************80
%
%% CH_EQI is a case insensitive comparison of two characters for equality.
%
%  Example:
%
%    CH_EQI ( 'A', 'a' ) is TRUE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 July 2000
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character C1, C2, the characters to compare.
%
%    Output, logical TRUEFALSE, is TRUE (1) if the characters are equal.
%
  if ( ch_cap ( c1 ) == ch_cap ( c2 ) )
    truefalse = 1;
  else
    truefalse = 0;
  end

  return
end
